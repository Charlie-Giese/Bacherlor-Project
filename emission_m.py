

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

parser = argparse.ArgumentParser(description='Plot an image in fits format.')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")
parser.add_argument("-c", "--contour_bias", help="value of contour level bias", default=False)
#parser.add_argument("-d", "--d_range", nargs="+", help="range of data to plot", default=False)
parser.add_argument("-o", "--outputfile", help="name of output .png file (optional)", default=False)
#parser.add_argument("-e", "--emission_measure", help="Optional emission measure plot, True/False", default=False)
args = parser.parse_args()

inputfile=args.inputfile
contour_bias=args.contour_bias
#d_range=args.d_range
outputfile=args.outputfile
#emission_measure=args.emission_measure
 


def import_fits(filename):    #, d_range):
	"""Imports a FITS radio image from CASA pipeline"""
	
	fits_image_filename = os.getcwd()+filename
	with fits.open(fits_image_filename) as hdul:
		#print(hdul.info())
		data=hdul[0].data

		
	
	#now convert to Jy/str surface brigthness
	beam_fwhm = 7.2*u.arcsec
	fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
	beam_sigma = beam_fwhm * fwhm_to_sigma
	omega_B = 2 * np.pi * beam_sigma**2
	SB=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(omega_B)))[0,0,:,:].value
	
	"""SB is a 2-d array with Surface_brigthness values in MJy/sr. For masked arrays, some values are nan, need to create a mask to 
	these values. Before that, need to set all values not within data range to NAN. Using SB.data can access the original array."""
	
	
	
	sb_maskedNAN = np.ma.array(SB, mask=np.isnan(SB))	
	#if d_range == True:
	#	sb_maskedNAN[(sb_maskedNAN < np.float32(d_range[0])) & (sb_maskedNAN > np.float32(d_range[1]))] = np.nan
	
	array = sb_maskedNAN
	
	return array
		

def em(inputfile, contour_bias, outputfile):  #,d_range)
	S=import_fits(inputfile)  #, d_range)

	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = S * 10**-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = -1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1	
	
	c=contour_bias
	
	hdu = fits.open(os.getcwd()+inputfile)
	hrd=hdu[0].header
	wcs = WCS(hrd)
	

	norm = ImageNormalize(emission_measure, interval=MinMaxInterval(),
                      stretch=SqrtStretch())
	
	fig=plt.figure(1) 
	ax=fig.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	em_map=ax.imshow(emission_measure, origin='lower', cmap='plasma')
	ax.set_xlabel('Right Ascension')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Emission Measure')
	if contour_bias != False:
		n = 7 #number of contours
		ub = np.max(emission_measure) #- (p-1)/(n)
		lb = np.min(emission_measure) #+ (p-1)/(n)
		spacing = c
		def level_func(lb, ub, n, spacing=1.1):
			span = (ub-lb)
			dx = 1.0 / (n-1)
			return [lb + (i*dx)**spacing*span for i in range(n)]
		levels=level_func(lb, ub, n, spacing=float(c))
		print(levels)
		ax.contour(emission_measure, levels=levels, colors='white', alpha=0.5)
	plt.show()
	
result=em(inputfile, contour_bias, outputfile)















