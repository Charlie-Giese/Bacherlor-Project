#Author: Charlie Giese, Created: 17/2/21

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
import astropy.units as u
import math as m
import os
from astropy.coordinates import SkyCoord
import sys, getopt
from astropy.wcs import WCS
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
from matplotlib import cm
import argparse
from radio_beam import Beam
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize)
from astropy.visualization.wcsaxes import Quadrangle
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.nddata import Cutout2D
from scipy import ndimage, misc

fontsize=14
font = {'family' : 'DejaVu Sans',
'size' : fontsize}

class MidpointNormalize(mpl.colors.Normalize):
    """
    class to help renormalize the color scale
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))



images = [
'smoothed_fits/standard/unmasked/4-5_smoothed.fits',
'smoothed_fits/standard/unmasked/5-6_smoothed.fits',
'smoothed_fits/standard/unmasked/6-7_smoothed.fits',
'smoothed_fits/standard/unmasked/7-8_smoothed.fits',
'smoothed_fits/standard/unmasked/8-9_smoothed.fits',
'smoothed_fits/standard/unmasked/9-10_smoothed.fits',
'smoothed_fits/standard/unmasked/10-11_smoothed.fits',

]

rms = [
0.0001, # 4-5
0.00015, # 5-6
0.00013, # 6-7
0.0001, # 7-8
0.0001, # 8-9
0.0001, # 9-10
0.0001, # 10-11
]


def s_nu(nu,s0,alpha): # f(xdata,a0,a1)
  return s0*nu**alpha

arrays=[]
i=0

freq_ar = np.array([4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5])

for im in images:
	with fits.open(im) as hdul:
		data=hdul[0].data[0,0,:,:]
	header = fits.getheader(im)
	wcs=WCS(header)
	print(wcs)
	beam = Beam.from_fits_header(header)
	SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
	where_are_NaNs = np.isnan(SB)
	SB[where_are_NaNs] = 0.0
	#SB_masked = np.ma.masked_invalid(SB)
	where_are_neg = np.where(SB < 0)
	SB[where_are_neg] = 0.0
	result = ndimage.zoom(SB, 1/4.0)
	print(np.shape(result))
	i+=1
	arrays.append(result)

spectral_ar = np.ones_like(arrays[0])*10

for freq in freq_ar:
	for i in range(np.shape(spectral_ar)[0]-1):
		for j in range(np.shape(spectral_ar)[1]-1):
			fluxes = []


			fluxes.append(arrays[0][i,j])
			fluxes.append(arrays[1][i,j])
			fluxes.append(arrays[2][i,j])
			fluxes.append(arrays[3][i,j])
			fluxes.append(arrays[4][i,j])
			fluxes.append(arrays[5][i,j])
			fluxes.append(arrays[6][i,j])

			flux_ar = np.array(fluxes)
			indices = np.where(flux_ar > 0.)
			flux_final=flux_ar[indices]
			freq_final = freq_ar[indices]

			if np.size(flux_final) >= 3:
				yerr = rms
				try:
					coef, cov = curve_fit(s_nu,freq_final,flux_final)#,sigma=yerr,absolute_sigma=True)
					s0        = coef[0]
					#freq_err  = np.sqrt(cov[0,0]) # [[sxx sxy][syx syy]]
					alpha     = coef[1]
					print(alpha)
					#alpha_err = np.sqrt(cov[1,1])
					spectral_ar[i,j] = alpha
				except:
					continue

wcs_input_dict = {
    'CTYPE1': 'RA---SIN',
    'CUNIT1': 'deg',
    'CDELT1': -0.0003333333333333,
    'CRPIX1': 321.,
    'CRVAL1': 350.20125,
    'NAXIS1': 160,
    'CTYPE2': 'DEC--SIN',
    'CUNIT2': 'deg',
    'CDELT2': 0.0003333333333333,
    'CRPIX2': 321.,
    'CRVAL2': 61.20166666429,
    'NAXIS2': 160
}
wcs2 = WCS(wcs_input_dict)

where_out_bounds1 = np.where(spectral_ar > 1.)
spectral_ar[where_out_bounds1] = np.nan
where_out_bounds2 = np.where(spectral_ar < -1.)
spectral_ar[where_out_bounds2] = np.nan

norm = MidpointNormalize(midpoint = -0.5)
figure=plt.figure(num=1)

ax=figure.add_subplot(111, projection=wcs2, slices=('x','y'))

main_image=ax.imshow(X=spectral_ar, cmap='seismic', origin='lower', norm=norm)
cbar=figure.colorbar(main_image)
ax.set_xlabel('Right Ascension J2000', fontdict=font)
ax.set_ylabel('Declination J2000', fontdict=font)
cbar.set_label('\u03B1', fontdict=font)
ra = ax.coords[0]
ra.set_format_unit('degree', decimal=True)

dec=ax.coords[1]
dec.set_format_unit('degree', decimal=True)
plt.show()
