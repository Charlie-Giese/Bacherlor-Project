

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
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize)
import argparse
from radio_beam import Beam
from astropy.visualization.wcsaxes import SphericalCircle

parser = argparse.ArgumentParser(description='Plot an image in fits format.')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")
parser.add_argument("-d", "--d_range", nargs="+", help="range of data to plot", default=False)
parser.add_argument("-s", '--imsize', help='Size of cutout image_plot', default=400)
parser.add_argument("-o", "--outputfile", help="name of output .png file (optional)", default=False)
args = parser.parse_args()
inputfile=args.inputfile

d_range=args.d_range
imsize=float(args.imsize)
outputfile=args.outputfile

#configuring font
fontsize=11
font = {'family' : 'DejaVu Sans',
'size' : fontsize}

### IMAGE PLOT CODE; INCLUDES A FUNCTION FOR IMPORTING FITS FILE, CONVERTING TO APPROPRIATE UNITS AND THEN PLOTTING ###

def jy_beam_MJy_sr(data, header):
	"""Funtion which converts fits data into units of Jy/sr"""

	if header['BUNIT'] == 'Jy/beam':
		beam = Beam.from_fits_header(header)
		SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
	elif header['BUNIT'] == 'JY/BEAM':
		beam = Beam.from_fits_header(header)
		SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
	elif header['BUNIT'] == 'MJy/sr':
		SB=data

	return SB


def import_fits(filename, d_range):
	"""Imports a FITS radio image from CASA pipeline and returns the image array in units of MJy/sr and the header"""

	fits_image_filename = filename
	with fits.open(fits_image_filename) as hdul:
		data=hdul[0].data
	header = fits.getheader(fits_image_filename)

	ndims = data.ndim
	if ndims == 4:
		image=data[0,0,:,:]
	else:
		image=data

	ar_u = jy_beam_MJy_sr(image, header)
	ar_c=mask(ar_u)
	if d_range != False:
		array=clip(ar_c, d_range)
		ar_c=mask(array)


	return ar_c,header


def mask(array):
	"""Masks the array so invalid values are discarded"""

	unmasked=np.array(array, dtype=np.float32)
	masked=np.ma.masked_invalid(unmasked)
	output=masked
	return output


def clip(array, d_range):
	"""Clips values outside of the data range, replaces them with max and min values"""


	array[(array < np.float32(d_range[0]))] = d_range[0]
	array[(array > np.float32(d_range[1]))] = d_range[1]

	array=mask(array)

	return array


def image_plot(inputfile, d_range, imsize, outputfile):
	"""Takes a scaled image in the form of a numpy array and plots it along with coordinates and a colorbar"""

	block = import_fits(inputfile, d_range)
	data = block[0]
	hrd = block[1]
	wcs = WCS(hrd)
	if hrd['TELESCOP'] != 'Spitzer':
		beam = Beam.from_fits_header(hrd)

	print(np.shape(data))
	norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=SqrtStretch())

	figure=plt.figure(num=1)
	if wcs.pixel_n_dim == 4:
		ax=figure.add_subplot(111, projection=wcs, slices=('x','y', 0, 0))
	elif wcs.pixel_n_dim ==2:
		ax=figure.add_subplot(111, projection=wcs, slices=('x','y'))


	main_image=ax.imshow(X=data, cmap='plasma', origin='lower', norm=norm, vmax=np.max(data) , vmin=np.min(data))
	cbar=figure.colorbar(main_image)

	if hrd['TELESCOP'] == 'Spitzer':
		ax.invert_xaxis()
		ax.invert_yaxis()

	#position=SkyCoord('23h20m48s', '+61d12m06s', frame='icrs')
	#centre=wcs.world_to_pixel(position)
	dims=np.shape(data)
	centre=(dims[0]/2., dims[1]/2.)
	ax.set_xlim(centre[0]-imsize, centre[0]+imsize)
	ax.set_ylim(centre[1]-imsize, centre[1]+imsize)



	ax.set_xlabel('Right Ascension J2000', fontdict=font)
	ax.set_ylabel('Declination J2000', fontdict=font)
	cbar.set_label('Surface Brigthness (MJy/Sr)', fontdict=font)
	ra = ax.coords[0]
	ra.set_format_unit('degree', decimal=True)

	dec=ax.coords[1]
	dec.set_format_unit('degree', decimal=True)

	#ax.set_title('Spitzer 24\u03BCm', fontdict=font)
	ax.set_title('8-12 GHz, 5\u03C3 mask', fontdict=font)

	if hrd['TELESCOP'] != 'Spitzer':
		beam = Beam.from_fits_header(hrd)
		c = SphericalCircle((350.34, 61.13)*u.degree, beam.major, edgecolor='black', facecolor='none',
	           				transform=ax.get_transform('fk5'))
		ax.add_patch(c)

	if outputfile != False:
		plt.savefig(os.getcwd()+'/thesis_figs/'+outputfile)
	plt.show()


result=image_plot(inputfile, d_range, imsize,  outputfile)
