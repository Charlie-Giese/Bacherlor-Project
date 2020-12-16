

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

parser = argparse.ArgumentParser(description='Plot an image in fits format.')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")
parser.add_argument("-cf", "--contour_file", help="fits file to plot contours from", default=False)
parser.add_argument("-ct", "--contour_type", help="automatic, sigma or manual") 
parser.add_argument("-cb", "--contour_bias", help="Provide value of contour level bias for automatic contours")
parser.add_argument("-cl", "--contour_levels", nargs="+", help="Provide sigma levels to plot contours at for sigma contours or MJy/sr values for manual contours")
parser.add_argument("-d", "--d_range", nargs="+", help="range of data to plot", default=False)
parser.add_argument("-o", "--outputfile", help="name of output .png file (optional)", default=False)
args = parser.parse_args()
inputfile=args.inputfile

contour_file=args.contour_file
contour_type=args.contour_type
contour_bias=args.contour_bias

d_range=args.d_range
outputfile=args.outputfile



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
		
	
def image_plot(inputfile, d_range, outputfile):
	"""Takes a scaled image in the form of a numpy array and plots it along with coordinates and a colorbar"""
	
	block = import_fits(inputfile, d_range)
	data = block[0]
	hrd = block[1]
	wcs = WCS(hrd)
	if hrd['TELESCOP'] != 'Spitzer':
		beam = Beam.from_fits_header(hrd)
		

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
		

	coord = SkyCoord('23h23m37s +61:04, frame='J200')
		
		
	ax.set_xlabel('Right Ascension J200')
	ax.set_ylabel('Declination J2000')
	cbar.set_label('Surface Brigthness (MJy/Sr)')
	
	if contour_file != False:
		contours=contour_plot(ax, contour_file, contour_type, beam)
		
		
	plt.show()
	if outputfile != False:
		plt.savefig(os.getcwd()+outputfile)


### CODE FOR IMPORTING A FITS FILE FOR CONTOURS AND THEN OVERPLOTTING THOSE CONTOURS ON THE IMAGE CREATED ABOVE ###



def import_contours(filename):
	"""Imports a FITS radio image from CASA pipeline, and coverts units to MJy/sr if not already in those units."""
	
	contour_fits = filename
	with fits.open(contour_fits) as hdul:
		data=hdul[0].data
	header = fits.getheader(contour_fits)
	
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



def contour_plot(ax, contour_file, contour_type, image_beam):
	"""Plots contours over the main image, can be either from the same image or from a different one, fits file."""

	cblock=import_contours(contour_file)
	cdata=cblock[0]
	cheader=cblock[1]
	cwcs = WCS(cheader)

	contour_type=str(contour_type)
	ax.set_autoscale_on(False)

	norm = ImageNormalize(cdata, interval=MinMaxInterval(), stretch=SqrtStretch())

	if contour_type == "automatic":
		spacing = contour_bias
		n = 7 #number of contours
		ub = np.max(contour_data) 
		lb = np.min(contour_data) 
		def level_func(lb, ub, n, spacing=1.1):
			span = (ub-lb)
			dx = 1.0 / (n-1)
			return [lb + (i*dx)**spacing*span for i in range(n)]
		levels=level_func(lb, ub, n, spacing=float(spacing))
		print('Generated Levels for contour are: ', levels, 'MJy/sr')
		ax.contour(cdata, levels=levels, colors='white', transform=ax.get_transform(cwcs.celestial), alpha=0.5)
	if contour_type == "sigma":
		contour_levels=list(map(float, args.contour_levels))
		sigma_levels=np.array(contour_levels)
		rms=(0.00004*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(image_beam))
		levels=rms*sigma_levels
		print('Sigma Levels for contour are: ', sigma_levels, 'sigma')
		ax.contour(cdata, levels=levels, norm=norm, transform=ax.get_transform(cwcs.celestial), colors='white', alpha=0.5)
	if contour_type == "manual":
		levels=list(map(float, args.contour_levels))
		print('Contour Levels are: ', levels, 'MJy/sr')
		ax.contour(cdata, levels=levels, norm=norm, transform=ax.get_transform(cwcs.celestial), colors='white', alpha=0.5)


	

result=image_plot(inputfile, d_range,  outputfile)









