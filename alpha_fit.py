#Author: Charlie Giese, Created: 17/2/21

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
from astropy.coordinates import SkyCoord
import sys, getopt
from astropy.wcs import WCS
from scipy.optimize import curve_fit
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

coord_ll = SkyCoord(350.30, 61.15, unit='deg', frame='fk5')
coord_ur = SkyCoord(350.10, 61.25, unit='deg', frame='fk5')

arrays=[]
i=0

freqs = [4.486956471109,
 		 5.510946533326,
		 6.4869370613760005,
		 7.510927123593,
		 8.486918282031999,
		 9.510908420254,
		 10.48689902075]

for im in images:
	with fits.open(im) as hdul:
		data=hdul[0].data[0,0,:,:]
	header = fits.getheader(im)
	wcs=WCS(header)[0,0,:,:]
	index_ll = wcs.world_to_array_index(coord_ll)
	index_ur = wcs.world_to_array_index(coord_ur)

	beam = Beam.from_fits_header(header)
	SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value

	ar = SB[index_ll[0] : index_ur[0], index_ll[1] : index_ur[1]]

	arrays.append(ndimage.zoom(ar, 4.0))
	i+=1


(330, 289)

spectral_ar = np.empty_like(arrays[0])

for freq in freqs:
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

			x = np.sum(np.array(fluxes))
			if np.isnan(x) == False:
				yerr = rms
				coef, cov = curve_fit(s_nu,freqs,fluxes,sigma=yerr,absolute_sigma=True)
				s0        = coef[0]
				freq_err  = np.sqrt(cov[0,0]) # [[sxx sxy][syx syy]]
				alpha     = coef[1]
				print(alpha)
				alpha_err = np.sqrt(cov[1,1])
				spectral_ar[i,j] = alpha





norm = ImageNormalize(spectral_ar, interval=MinMaxInterval(), stretch=SqrtStretch())
figure=plt.figure(num=1)

ax=figure.add_subplot(111, projection=wcs, slices=('x','y'))

main_image=ax.imshow(X=spectral_ar, cmap='plasma', origin='lower', norm=norm)
cbar=figure.colorbar(main_image)
ax.set_xlabel('Right Ascension J2000', fontdict=font)
ax.set_ylabel('Declination J2000', fontdict=font)
cbar.set_label('Surface Brigthness (MJy/Sr)', fontdict=font)
ra = ax.coords[0]
ra.set_format_unit('degree', decimal=True)

dec=ax.coords[1]
dec.set_format_unit('degree', decimal=True)
plt.show()
