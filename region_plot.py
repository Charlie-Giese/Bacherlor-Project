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
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.nddata import Cutout2D
import matplotlib


#configuring matplotlib
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['figure.figsize'] = [6.8,5.5]
matplotlib.rcParams['figure.dpi'] = 120
matplotlib.rcParams['font.sans-serif'] = "Nimbus Roman"

#coordinates of BD+60 2522
star_coord = SkyCoord("23h20m44.5s +61d11m40.5s", frame = 'icrs')


#fits_image_filename = 'fits/ngc7635_g0.5_briggs1.0_nsig5.image.5sig_masked.fits'
fits_image_filename = 'fits/ngc7635_g0.5_briggs1.0_nsig5.image.tt0.fits'

with fits.open(fits_image_filename) as hdul:
	data=hdul[0].data[0,0,:,:]
header = fits.getheader(fits_image_filename)

beam = Beam.from_fits_header(header)
SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
SB_masked=np.ma.masked_invalid(SB)

def optical_depth(flux):

	SB_erg=flux * 1e-17
	T=8000 #K
	nu=8e9 # central frequency of our broadband observations
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	tau = - np.log( 1 - (c**2 *SB_erg)/(2*k_b*T*(nu**2)))
	return tau



# Position of NGC 7635: 23 20 48.3 +61 12 06
#position = SkyCoord('23h20m48s', '+61d12m06s', frame='fk5')
# subset of image, centred on NGC7635
#cutout = Cutout2D(SB_masked, position, (1000, 1000), wcs=wcs)

plot_data = SB_masked

wcs = WCS(header, naxis=2)
norm = ImageNormalize(plot_data, interval=MinMaxInterval(), stretch=SqrtStretch())

figure=plt.figure(num=1)
ax=figure.add_subplot(111, projection=wcs)
main_image=ax.imshow(X=plot_data, cmap='plasma', origin='lower', norm=norm)	#, vmax=np.max(data) , vmin=np.min(data))
cbar=figure.colorbar(main_image)

star_index = wcs.world_to_pixel(star_coord)
star=ax.scatter(star_index[0], star_index[1], marker='*', c='w')

ax.set_xlabel('Right Ascension J2000')
ax.set_ylabel('Declination J2000')
cbar.set_label('Surface Brigthness (MJy/Sr)')

dims=np.shape(plot_data)
centre=(dims[0]/2., dims[1]/2.)
ax.set_xlim(centre[0]-300, centre[0]+300)
ax.set_ylim(centre[1]-300, centre[1]+300)


#regions = open('tau_region_deg.txt', 'r')
regions = open('regions_degrees.txt', 'r')
#350.199, 61.2321, 0.0230107, 0.0246958
i=0
labels=['a','b','c','d']

for reg in regions.readlines():
	region = reg.split(",")
	h = float(region[3])
	y = (float(region[1])-0.5*h)
	w = (float(region[2])/np.cos(y*np.pi/180.))
	x = (float(region[0])-0.5*w)

	r = Quadrangle	((x, y)*u.degree, w*u.degree, h*u.degree, vertex_unit='degree',
              label='labels[i]', edgecolor='white', facecolor='none', linestyle='-',
              transform=ax.get_transform('fk5'))


	ax.add_patch(r)

	plt.text(float(region[0]), float(region[1]), labels[i], transform=ax.get_transform('fk5'), c='w')
	i+=1

ra = ax.coords[0]
ra.set_format_unit('degree', decimal=True)

dec=ax.coords[1]
dec.set_format_unit('degree', decimal=True)

dims=np.shape(SB_masked)
centre=(dims[0]/2, dims[1]/2)
size=400
ax.set_xlim(centre[0]-size, centre[0]+size)
ax.set_ylim(centre[1]-size, centre[1]+size)

beam = Beam.from_fits_header(header)
c = SphericalCircle((350.34, 61.13)*u.degree, beam.major, edgecolor='white', facecolor='none',
           transform=ax.get_transform('fk5'))
ax.add_patch(c)

plt.show()
