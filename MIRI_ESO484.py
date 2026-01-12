"""
Stage 1, 2, and 3 MIRI pipeline reduction for The Cosmic Eye
https://www.stsci.edu/jwst-program-info/download/jwst/pdf/4125/

"""

# _______________________________________________________________________________________________#
import os
import pathlib
import shutil

# Set up CRDS path and server environment variables
# Note: the CRDS_PATH and CRDS_SERVER_URL *MUST* be set before importing the jwst module.
os.environ["CRDS_PATH"] = "calibration_reference_files"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"

import stcal
import jwst
from jwst.pipeline import Detector1Pipeline  # stage 1
from jwst.pipeline.calwebb_image2 import Image2Pipeline  # stage 2
from jwst.pipeline.calwebb_image3 import Image3Pipeline  # stage 3

from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

# _______________________________________________________________________________________________#
#                                        Set up run Parameters
# _______________________________________________________________________________________________#

target = 'ESO484'
obs_list = ['02101_1', '02101_3', '02103_1', '02103_3', '02105_1', '02105_3']

file_dir = pathlib.Path(__file__).parent.resolve()

# For future reference, keep track of the current JWST and STCAL versions
print(f"# JWST pipe version = {jwst.__version__}")
print(f"# STCAL version = {stcal.__version__}")

for obs in obs_list:
    # raw_dir is where the uncalibrated (*uncal.fits) are, which were downloaded from MAST
    raw_dir = file_dir / f"{target}/MIRI/{obs}/raw"
    reduced_dir = file_dir / f"{target}/MIRI/{obs}/reduced"

    reduce_stage1 = False
    reduce_stage2 = True
    reduce_stage3 = True

    create_config_file = False  # Turn this off after fist creation if want to change defaults
    # _______________________________________________________________________________________________#
    #                                        Directory Setup
    # _______________________________________________________________________________________________#

    # Cosmic Eye only only has sci files without separate sky exposures.

    # set up reduction directories
    stage1_dir = reduced_dir / "stage1"
    stage2_dir = reduced_dir / "stage2"
    stage3_dir = reduced_dir / "stage3"

    # if the directories don't already exist, make them
    for folder in [raw_dir, stage1_dir, stage2_dir, stage3_dir]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    # Preamble: copy params files
    for stage in [1, 2, 3]:
        src = file_dir / f"{target}/MIRI/stage{stage}_params.asdf"
        dst_dir = eval(f"stage{stage}_dir")
        dst = dst_dir / f"stage{stage}_params.asdf"
        if not dst.exists():
            shutil.copy(src, dst)

    # ------------------------------ Create the config files ----------------------------------------#
    # The default values in these files can be overwritten there, or at the call to the pipeline by
    # passing in python dictionaries using the keywords in the config files.
    if create_config_file:
        pipelines = [Detector1Pipeline(), Image2Pipeline(), Image3Pipeline()]
        directories = [stage1_dir, stage2_dir, stage3_dir]
        for i, (pipeline, directory) in enumerate(zip(pipelines, directories)):
            config_filename = directory / f"stage{i + 1}_params.asdf"
            pipeline.export_config(config_filename)

    # _______________________________________________________________________________________________#
    #                                         Stage 1
    # _______________________________________________________________________________________________#

    if reduce_stage1 is True:

        # get the files to work on.
        raw_uncal_files = sorted(raw_dir.glob("*_uncal.fits"))
        print(f"******** Step 1: Working on {len(raw_uncal_files)} uncal files: ******** \n")
        print(raw_uncal_files)

        # run stage 1 of the pipeline
        for uncal_file in raw_uncal_files:
            print("Applying Stage 1 Corrections & Calibrations to: " + uncal_file.name)

            result = Detector1Pipeline.call(
                uncal_file,
                save_results=True,
                output_dir=str(stage1_dir),
                config_file=str(stage1_dir / "stage1_params.asdf"),
            )

    # _______________________________________________________________________________________________#
    #                                           Stage 2
    # _______________________________________________________________________________________________#

    if reduce_stage2 is True:

        # get the intermediate rate files
        inter_rate_files = sorted(stage1_dir.glob("*_rate.fits"))
        print(f"******** Step 2: Working on {len(inter_rate_files)} rate files: ******** \n")
        print(inter_rate_files)

        # Run the stage 2 of the pipeline
        for rate_file in inter_rate_files:
            print("Applying Stage 2 Calibrations & Corrections to: " + rate_file.name)

            result = Image2Pipeline.call(
                rate_file,
                save_results=True,
                output_dir=str(stage2_dir),
                config_file=str(stage2_dir / "stage2_params.asdf"),
            )

    # _______________________________________________________________________________________________#
    #                                        Stage 3
    # _______________________________________________________________________________________________#

    if reduce_stage3 is True:

        # get the intermediate cal files
        inter_cal_files = [str(filename) for filename in sorted(stage2_dir.glob("*_cal.fits"))]
        print(f"******** Step 3: Working on {len(inter_cal_files)} cal files: ******** \n")
        print(inter_cal_files)

        # creat the JSON file describing the file associations
        out_asn = asn_from_list.asn_from_list(
            items=inter_cal_files, rule=DMS_Level3_Base, product_name="stage3"
        )
        output_asn = stage3_dir / f"{target}_calwebb3.json"

        with open(output_asn, "w") as outfile:
            name, serialized = out_asn.dump(format="json")
            outfile.write(serialized)

        # run the final stage
        result = Image3Pipeline.call(
            str(output_asn),  # Association (ASN) file listing the input exposures
            save_results=True,  # Write outputs of each step to disk
            output_dir=str(stage3_dir),  # Directory where outputs will be saved
            config_file=str(stage3_dir / "stage3_params.asdf"),
        )