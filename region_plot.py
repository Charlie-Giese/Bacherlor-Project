import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
from astropy.coordinates import SkyCoord
import sys, getopt
from astropy.wcs import WCS
from matplotlib import cm
import argparse
from radio_beam import Beam
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize)
from astropy.visualization.wcsaxes import Quadrangle
from matplotlib.patches import Rectangle




fits_image_filename = 'fits/ngc7635_g0.5_briggs1.0_nsig5.image.tt0.fits'

with fits.open(fits_image_filename) as hdul:
	data=hdul[0].data[0,0,:,:]
header = fits.getheader(fits_image_filename)

beam = Beam.from_fits_header(header)
SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
SB_masked=np.ma.masked_invalid(data)

wcs = WCS(header)
norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=SqrtStretch())

figure=plt.figure(num=1)
ax=figure.add_subplot(111, projection=wcs, slices=('x','y', 0, 0))
main_image=ax.imshow(X=SB_masked, cmap='plasma', origin='lower', norm=norm)#, vmax=np.max(data) , vmin=np.min(data))
cbar=figure.colorbar(main_image)

regions = open('regions_degrees.txt', 'r')

i=0
labels=['a','b','c','d','e']

for reg in regions.readlines():
	region = reg.split(",")
	x = (float(region[0]) - float(region[2])/2)*u.degree
	y = (float(region[1]) - float(region[3])/2)*u.degree
	w = float(region[3])*u.degree
	h = float(region[2])*u.degree


	r = Quadrangle	((x, y), w, h,
              label='labels[i]', edgecolor='white', facecolor='none', linestyle='-',
              transform=ax.get_transform('fk5'))
	ax.add_patch(r)

	plt.text(float(region[0]), float(region[1]), labels[i], transform=ax.get_transform('fk5'), c='w')
	i+=1

ax.set_xlabel('Right Ascension J200')
ax.set_ylabel('Declination J2000')
cbar.set_label('Surface Brigthness (MJy/Sr)')


plt.show()
