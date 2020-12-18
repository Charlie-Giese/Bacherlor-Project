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
from matplotlib.patches import Rectangle

fits_image_filename = 'fits/ngc7635_Xband_wide_g0.5_briggs1.0_1000iter_interactive.image.tt0.fits'

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

r = Rectangle((355*u.degree, +61.12*u.degree), 1.397*u.arcmin , 0.63*u.arcmin,
              label='Rectangle', edgecolor='red', facecolor='none', linestyle='--',
              transform=ax.get_transform('icrs'))
ax.add_patch(r)



ax.set_xlabel('Right Ascension J200')
ax.set_ylabel('Declination J2000')
cbar.set_label('Surface Brigthness (MJy/Sr)')

plt.show()
