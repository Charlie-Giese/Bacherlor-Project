

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
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
parser.add_argument("-ct", "--contour_type", help="automatic (contour_bias) or defined (contour_levels)") 
parser.add_argument("-cb", "--contour_bias", help="value of contour level bias")
parser.add_argument("-cl", "--contour_levels", nargs="+", help="sigma levels to plot contours at")
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


def import_fits(filename, d_range):
	"""Imports a FITS radio image from CASA pipeline"""
	
	fits_image_filename = os.getcwd()+filename
	with fits.open(fits_image_filename) as hdul:
		#print(hdul.info())
		data=hdul[0].data

	header = fits.getheader(fits_image_filename)
	if header['BUNIT'] == 'Jy/beam':
		beam = Beam.from_fits_header(header)
		SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam)))[0,0,:,:].value
	elif header['BUNIT'] == 'JY/BEAM':
		beam = Beam.from_fits_header(header)
		SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam)))[0,0,:,:].value
	elif header['BUNIT'] == 'MJy/sr':
		SB=data

	
	"""SB is a 2-d array with Surface_brigthness values in MJy/sr. For masked arrays, some values are nan, need to create a mask to 
	these values. Before that, need to set all values not within data range to NAN. Using SB.data can access the original array."""
	
	sb_maskedNAN = np.ma.array(SB, mask=np.isnan(SB))	
	if d_range != False:
		sb_maskedNAN[(sb_maskedNAN < np.float32(d_range[0]))] = d_range[0]
		sb_maskedNAN[(sb_maskedNAN > np.float32(d_range[1]))] = d_range[1]
	
	array = np.ma.array(sb_maskedNAN, mask=np.isnan(sb_maskedNAN))	
	
	return array
	
def image_plot(inputfile, d_range, outputfile):
	"""Takes a scaled image in the form of a numpy array and plots it along with coordinates and a colorbar"""
	
	hdu = fits.open(os.getcwd()+inputfile)
	hrd=hdu[0].header
	wcs = WCS(hrd)
	
	data=import_fits(inputfile, d_range)
	
	#if d_range != False:
	norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=SqrtStretch())
	#else:
	#	norm = ImageNormalize(data, vmin=d_range[0], vmax=d_range[1], stretch=SqrtStretch())
	
	figure=plt.figure(num=1)
	if wcs.pixel_n_dim == 4:
		ax=figure.add_subplot(111, projection=wcs, slices=('x','y', 0, 0))
	elif wcs.pixel_n_dim ==2:
		ax=figure.add_subplot(111, projection=wcs, slices=('x','y'))
	main_image=ax.imshow(X=data, cmap='plasma', origin='lower', norm=norm, vmax=np.max(data)- 5 , vmin=np.min(data))
	cbar=figure.colorbar(main_image)
	if hrd['TELESCOP'] == 'Spitzer':
		ax.invert_xaxis()
		ax.invert_yaxis()
	ax.set_xlabel('Right Ascension J2000')
	ax.set_ylabel('Declination J2000')
	cbar.set_label('Surface Brigthness (MJy/Sr)')
	if contour_file != False:
		image_beam = Beam.from_fits_header(hrd)
		contours=contour_plot(ax, contour_file, contour_type, image_beam)
	plt.show()
	if outputfile != False:
		plt.savefig(os.getcwd()+outputfile)


### CODE FOR IMPORTING A FITS FILE FOR CONTOURS AND THEN OVERPLOTTING THOSE CONTOURS ON THE IMAGE CREATED ABOVE ###



def import_contours(filename):
	"""Imports a FITS radio image from CASA pipeline, and coverts units to MJy/sr if not already in those units."""
	
	fits_image_filename = os.getcwd()+filename
	with fits.open(fits_image_filename) as hdul:
		#print(hdul.info())
		data=hdul[0].data

	header = fits.getheader(fits_image_filename)
	if header['BUNIT'] == 'Jy/beam':
		beam = Beam.from_fits_header(header)
		SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam)))[0,0,:,:].value
	elif header['BUNIT'] == 'JY/BEAM':
		beam = Beam.from_fits_header(header)
		SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam)))[0,0,:,:].value
	elif header['BUNIT'] == 'MJy/sr':
		SB=data*(u.MJy/u.sr)
	"""SB is a 2-d array with Surface_brigthness values in MJy/sr. For masked arrays, some values are nan, need to create a mask to 
	these values. Before that, need to set all values not within data range to NAN. Using SB.data can access the original array."""
	print('The contour data is in units of ', SB.unit)
	SB=SB.value
	
	array = np.ma.array(SB, mask=np.isnan(SB))	
	
	return array



def contour_plot(ax, contour_file, contour_type, beam):
	"""Plots contours over the main image, can be either from the same image or from a different one, fits file."""


	hdu_c = fits.open(os.getcwd()+contour_file)
	hrd_c=hdu_c[0].header
	wcs_c = WCS(hrd_c)
	contour_data=import_contours(contour_file)	
	print(np.max(contour_data))
	contour_type=str(contour_type)
	ax.set_autoscale_on(False)
	data=import_fits(contour_file, d_range)
	
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
		ax.contour(contour_data, levels=levels, colors='white', transform=ax.get_transform(wcs_c.celestial), alpha=0.5)
	if contour_type == "defined":
		contour_levels=list(map(float, args.contour_levels))
		sigma_levels=np.array(contour_levels)
		rms=(0.00004*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))
		levels=rms*sigma_levels
		print('Sigma Levels for contour are: ', sigma_levels, 'sigma')
		ax.contour(contour_data, levels=levels, transform=ax.get_transform(wcs_c.celestial), colors='white', alpha=0.5)


	

result=image_plot(inputfile, d_range,  outputfile)
















