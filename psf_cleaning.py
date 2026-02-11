import sys
import stpsf
import matplotlib.pylab as plt
from astropy.nddata import Cutout2D
from reproject import reproject_interp
import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils.aperture import CircularAperture, ApertureStats
from astropy.wcs import WCS
import matplotlib.colors as colors
from scipy.stats import linregress
from skimage import restoration
from clij2fft.richardson_lucy import richardson_lucy_nc
import asdf

def main(instrument,filter_name):
    # Load parameters
    with asdf.open('psf_cleaning_params.asdf') as af:
        params = af['params']
    
    if filter_name not in params:
        print(f"Filter {filter_name} not found in parameters file.")
        sys.exit(1)
    
    p = params[filter_name]
    subarray = p['subarray']
    bg_coords = p['bg_coords']
    sn_threshold = p['SN']
    position = p['position']
    size_pixels = p['size_pixels']
    jitter_sigma = p['jitter_sigma']
    alpha = p['alpha']
    num_iter = p['num_iter']
    bg_coords_final = p['bg_coords_final']
    psf_fov = p['psf_fov']
        
    # File paths
    file_full = f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_FULL_stage3.fits'
    file_sub = f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_{subarray}_stage3.fits'
    
    # Load FULL data
    with fits.open(file_full) as hdul:
        data = hdul['SCI'].data
        header = hdul['SCI'].header
        weight_map = hdul['WHT'].data
    
    # Load SUB data
    with fits.open(file_sub) as hdul:
        data_replace = hdul['SCI'].data
        header_replace = hdul['SCI'].header
    
    # Reproject SUB to FULL
    wcs = WCS(header)
    wcs_replace = WCS(header_replace)
    data_replace_reproj, _ = reproject_interp((data_replace, wcs_replace), wcs)
    
    # Fit transfer function
    mask_sat_full = np.isnan(data)
    mask_noise_sub = (data_replace_reproj < 0.1)
    mask_anchor = (~mask_sat_full) & (~mask_noise_sub) & (~np.isnan(data_replace_reproj))
    
    pixels_full = data[mask_anchor]
    pixels_sub = data_replace_reproj[mask_anchor]
    
    slope, intercept, _, _, _ = linregress(pixels_sub, pixels_full)
    
    data_sub_calibrated = (data_replace_reproj * slope) + intercept
    
    # Merge
    data_merged = data.copy()
    replace_mask = np.isnan(data) & (~np.isnan(data_sub_calibrated))
    data_merged[replace_mask] = data_sub_calibrated[replace_mask]
    
    # Background subtraction
    mask = (weight_map == 0)
    apertures = CircularAperture(bg_coords, r=3)
    sigma_clip = SigmaClip(sigma=3.0)
    phot_stats = ApertureStats(data_merged, apertures, mask=mask, sigma_clip=sigma_clip)
    
    final_bg_value = np.nanmean(phot_stats.median)
    final_bg_noise = np.nanmean(phot_stats.std)
    
    data_subtracted = data_merged - final_bg_value
    
    # SN cut
    sn_map = data_subtracted / final_bg_noise
    data_clean_sn = np.where(sn_map < sn_threshold, np.nan, data_subtracted)
    
    # Cutout
    cutout = Cutout2D(data_merged, position, size_pixels, wcs=wcs)
    data_source = cutout.data
    cutout_wcs = cutout.wcs
    data_clean_for_image = Cutout2D(data_clean_sn, position, size_pixels).data
    
    # PSF
    inst = stpsf.setup_sim_to_match_file(file_full)
    inst.options['jitter_sigma'] = jitter_sigma
    print(inst.detector_position)
    # To avoid aliasing warning, calculate PSF with smaller FOV and pad to match image size
    single_stpsf = inst.calc_psf(fov_pixels=psf_fov)
    psf_kernel_temp = single_stpsf['DET_DIST'].data
    # Pad the PSF to match the image size
    pad_width = (size_pixels - psf_fov) // 2
    psf_kernel = np.pad(psf_kernel_temp, pad_width, mode='constant', constant_values=0)
    # Sharpen the PSF to increase contrast while preserving total flux
    psf_kernel = psf_kernel ** alpha
    psf_kernel = psf_kernel / np.sum(psf_kernel)  # Renormalize to preserve flux
    
    # -----------------------------------------------------------------------
    # 2. Prepare the Science Data (Safety Steps)
    # -----------------------------------------------------------------------
    # Richardson-Lucy (the deconvolution algorithm) strictly requires:
    #  1. No NaNs
    #  2. No Negative values (it assumes Poisson statistics)
    data_clean = data_source.copy()
    # A. Handle NaNs
    data_clean = np.nan_to_num(data_clean, nan=0.0)
    # B. Handle Negatives (Background noise)
    # Since your data is background-subtracted, you likely have negative noise pixels.
    # We add a "pedestal" (constant background) to make everything positive, 
    # deconvolve, and then remove it later.
    min_val = np.min(data_clean)
    pedestal = 0
    if min_val < 0:
        pedestal = np.abs(min_val) + 1.0  # Add a buffer to be safe
        print(f"Adding pedestal of {pedestal:.4f} to handle negative pixels.")
        data_clean += pedestal
    
    # Deconvolve
    deconvolved_image = restoration.richardson_lucy(
        data_clean, 
        psf_kernel, 
        num_iter=num_iter, 
        clip=False
    )

    # deconvolved_image = richardson_lucy_nc(
    #     data_clean,
    #     psf_kernel,
    #     numiterations=num_iter,
    #     regularizationfactor=0.0
    # )

    if pedestal > 0:
        deconvolved_image -= pedestal
    
    # Residual clean
    apertures_final = CircularAperture(bg_coords_final, r=3)
    phot_stats_final = ApertureStats(data_merged, apertures_final, mask=mask, sigma_clip=sigma_clip)
    final_bg_value_res = np.nanmean(phot_stats_final.median)
    final_bg_noise_res = np.nanmean(phot_stats_final.std)
    
    residual_subtracted = deconvolved_image - final_bg_value_res
    sn_map_res = residual_subtracted / final_bg_noise_res
    residual_clean = np.where(sn_map_res < sn_threshold, np.nan, residual_subtracted)
    
    # Save fits files
    new_header = header.copy()
    new_header.update(cutout_wcs.to_header())
    new_header['HISTORY'] = 'Pixel replacement: saturated FULL pixels replaced with calibrated SUB640 data'
    new_header['HISTORY'] = f'Background subtraction and S/N cleaning with threshold {sn_threshold}'
    new_header['HISTORY'] = f'Richardson-Lucy deconvolution with {num_iter} iterations'
    new_header['HISTORY'] = f'PSF kernel sharpened with alpha={alpha}'
    
    # data_clean_for_image
    hdu1 = fits.PrimaryHDU(data=data_clean_for_image, header=new_header)
    hdu1.header['COMMENT'] = 'Cleaned image: centered, bg subtracted, SN cut'
    hdu1.writeto(f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_data_clean.fits', overwrite=True)
    
    # deconvolved_image
    hdu2 = fits.PrimaryHDU(data=deconvolved_image, header=new_header)
    hdu2.header['COMMENT'] = 'Deconvolved image'
    hdu2.writeto(f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_deconvolved.fits', overwrite=True)
    
    # residual_clean
    hdu3 = fits.PrimaryHDU(data=residual_clean, header=new_header)
    hdu3.header['COMMENT'] = 'Deconvolved image with bg subtraction and SN cut'
    hdu3.writeto(f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_residual_clean.fits', overwrite=True)
    
    # Save PNG
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].imshow(data_clean_for_image, norm=colors.LogNorm(), origin='lower', cmap='coolwarm')
    ax[0].set_title(f'Original Data ({filter_name})')
    ax[1].imshow(residual_clean, norm=colors.LogNorm(), origin='lower', cmap='coolwarm')
    ax[1].set_title('Deconvolved Data')
    plt.savefig(f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_preview.png')
    plt.close()
    
    # Save PSF kernel image
    plt.figure(figsize=(6,6))
    plt.imshow(psf_kernel, norm=colors.LogNorm(), origin='lower', cmap='coolwarm')
    plt.title('PSF Kernel')
    plt.savefig(f'ESO484/{instrument}/finals/{filter_name}/{filter_name}_psf_kernel.png')
    plt.close()

params = {
    'F430M': {
        'subarray': 'SUB640',
        'bg_coords': [(847, 4363), (3727, 589), (313, 1075), (4531, 709)],
        'SN': 5.0,
        'position': (991, 1106),
        'size_pixels': 1200,
        'jitter_sigma': 0.00,
        'alpha': 1.0,
        'psf_fov': 800,
        'num_iter': 2,
        'bg_coords_final': [(171, 101), (681, 1009), (921,1011)]
    },
    # Add more filters as needed
    'F335M': {
        'subarray': 'SUB640',
        'bg_coords': [(847, 4363), (3727, 589), (313, 1075), (4531, 709)],
        'SN': 5.0,
        'position': (991, 1106),
        'size_pixels': 1200,
        'jitter_sigma': 0.01,
        'alpha': 1.0,
        'psf_fov': 800,
        'num_iter': 2,
        'bg_coords_final': [(171, 101), (681, 1009), (921,1011)]
    },
    'F150W': {
        'subarray': 'SUB640',
        'bg_coords': [(3595, 3076), (2269, 325), (4098, 1163), (644, 3247)],
        'SN': 7.0,
        'position': (2109, 2315),
        'size_pixels': 2000,
        'jitter_sigma': 0.01,
        'alpha': 1,
        'psf_fov': 1000,
        'num_iter': 2,
        'bg_coords_final': [(343, 817), (1483, 1195), (1144,109), (1534,1288)]
    },
    'F200W': {
        'subarray': 'SUB640',
        'bg_coords': [(3595, 3076), (2269, 325), (4098, 1163), (644, 3247)],
        'SN': 7.0,
        'position': (2109, 2315),
        'size_pixels': 2000,
        'jitter_sigma': 0.01,
        'alpha': 1,
        'psf_fov': 1200,
        'num_iter': 2,
        'bg_coords_final': [(343, 817), (1483, 1195), (1144,109), (1534,1288)]
    },
    'F770W': {
        'subarray': 'SUB256',
        'bg_coords': [(1001, 303), (847, 769), (400,910), (651, 1061)],
        'SN': 7.0,
        'position': (759, 590),
        'size_pixels': 800,
        'jitter_sigma': 0.01,
        'alpha': 3,
        'psf_fov': 800,
        'num_iter': 2,
        'bg_coords_final': [(750, 150), (200, 100), (750,700), (400,100)]
    },
    'F1130W': {
        'subarray': 'SUB256',
        'bg_coords': [(1001, 303), (847, 769), (400,910), (651, 1061)],
        'SN': 7.0,
        'position': (759, 590),
        'size_pixels': 800,
        'jitter_sigma': 0.01,
        'alpha': 3,
        'psf_fov': 800,
        'num_iter': 2,
        'bg_coords_final': [(750, 150), (200, 100), (750,700), (400,100)]
    },
    'F2100W': {
        'subarray': 'SUB256',
        'bg_coords': [(1001, 303), (847, 769), (400,910), (651, 1061)],
        'SN': 7.0,
        'position': (759, 590),
        'size_pixels': 800,
        'jitter_sigma': 0.02,
        'alpha': 5,
        'psf_fov': 800,
        'num_iter': 2,
        'bg_coords_final': [(541, 207), (661, 76), (645,51), (689,365)]
    }

}

tree = {'params': params}
af = asdf.AsdfFile(tree)
af.write_to('psf_cleaning_params.asdf')

if __name__ == '__main__':

    instrument = 'MIRI'  # Assume NIRCam for now
    filter_name = 'F2100W'  # Example filter

    main(instrument, filter_name)