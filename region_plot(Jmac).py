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
#from astropy.visualization.wcsaxes import Quadrangle
#from matplotlib.patches import Rectangle
from matplotlib.patches import Rectangle
from astropy.nddata import Cutout2D


fits_image_filename = 'SPITZER_M1_15052544_0000_8_E6328389_maic.fits'

hdu = fits.open(fits_image_filename)
fits.info(fits_image_filename)
data =fits.getdata(fits_image_filename,ext=0)
where_are_NaNs = np.isnan(data)
data[where_are_NaNs] = 0

#quit()
#  print(hdul[0])
#  data=hdul[0].data
header = fits.getheader(fits_image_filename)
#print(data)
#beam = Beam.from_fits_header(header)
#SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
#SB_masked=np.ma.masked_invalid(data)

wcs = WCS(header)
norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=SqrtStretch())

# Position of NGC 7635: 23 20 48.3 +61 12 06
position = SkyCoord('23h20m48s', '+61d12m06s', frame='icrs')
# subset of image, centred on NGC7635
cutout = Cutout2D(data, position, (250, 250), wcs=wcs)

figure=plt.figure(num=1)
# use wcs of cutout image
ax=figure.add_subplot(111, projection=cutout.wcs)
main_image=ax.imshow(cutout.data, cmap='plasma', origin='lower', norm=norm)
cbar=figure.colorbar(main_image)
ax.grid()
ax.set_xlabel('Right Ascension J2000')
ax.set_ylabel('Declination J2000')
cbar.set_label('Surface Brigthness (MJy/Sr)')


#ra = ax.coords[0]
#ra.set_format_unit('degree', decimal=True)
#dec=ax.coords[1]
#dec.set_format_unit('degree', decimal=True)
#ax.set_xlim(350.1,350.4)
#ax.set_ylim(61.2,61.3)

regions = open('regions_degrees.txt', 'r')

r=0
for reg in regions.readlines():
  region = reg.split(",")
#x = float(region[0])*u.degree
#y = float(region[1])*u.degree
#w = float(region[2])*u.degree
#h = float(region[3])*u.degree
# Rectangle can't handle units...
  h = float(region[3])
  y = float(region[1])-0.5*h
  w = float(region[2])/np.cos(y*np.pi/180.)
  x = float(region[0])-0.5*w

  print("Rectangle:",x,y,w,h)
  r = Rectangle	((x, y), w, h,
    label='Rectangle', edgecolor='white', facecolor='none', linestyle='-',
    transform=ax.get_transform('icrs'))

  ax.add_patch(r)


plt.show()
