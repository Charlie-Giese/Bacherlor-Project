

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

parser = argparse.ArgumentParser(description='Calculate Emission measure and Optical depth as function of position')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")

args = parser.parse_args()

inputfile=args.inputfile

 


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


def import_fits(filename):
	"""Imports a FITS radio image from CASA pipeline and returns the image array in units of MJy/sr and the header"""
	
	fits_image_filename = os.getcwd()+'/fits/'+filename
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
	
	return ar_c,header


def mask(array):
	"""Masks the array so invalid values are discarded"""
	
	unmasked=np.array(array, dtype=np.float32)
	masked=np.ma.masked_invalid(unmasked)
	output=masked
	return output

	
def optical_depth(inputfile):

	SB=import_fits(inputfile)[0]
	SB_erg=SB * 1e-17

	T=8000 #K
	nu=8e9 # central frequency of our broadband observations
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs


	tau = - np.log( 1 - (c**2 *SB_erg)/(2*k_b*T*(nu**2)))
	
	
	hrd = import_fits(inputfile)[1]
	wcs = WCS(hrd)
	

	norm = ImageNormalize(emission_measure, interval=MinMaxInterval(),
                      stretch=SqrtStretch())
	
	fig=plt.figure(1) 
	ax=fig.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	em_map=ax.imshow(tau, origin='lower', cmap='plasma', norm=norm)
	ax.set_xlabel('Right Ascension')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Optical Depth, kpc cm^-6')
	plt.show()
		

def em(inputfile):
	S=import_fits(inputfile)[0]

	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = S * 10**-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = (-1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1)/(1e3)
	
	
	
	hrd = import_fits(inputfile)[1]
	wcs = WCS(hrd)
	

	norm = ImageNormalize(emission_measure, interval=MinMaxInterval(),
                      stretch=SqrtStretch())
	
	fig=plt.figure(1) 
	ax=fig.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	em_map=ax.imshow(emission_measure, origin='lower', cmap='plasma', norm=norm)
	ax.set_xlabel('Right Ascension')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Emission Measure, kpc cm^-6')
	plt.show()
	
	

	
	
	
	
emission_measure=em(inputfile)
tau=optical_depth(inputfile)















