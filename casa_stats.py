import numpy as np
import argparse
import os
from astropy.io import fits
import astropy.units as u

parser = argparse.ArgumentParser(description='Returns min and max of image in MJy/sr')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")
args = parser.parse_args()
inputfile=args.inputfile
def import_fits(filename):
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
	
	
	array = np.ma.array(sb_maskedNAN, mask=np.isnan(sb_maskedNAN))	
	
	print('Max value in image is :', np.max(sb_maskedNAN), ' MJy/sr')
	print('Min value in image is :', np.min(sb_maskedNAN), ' MJy/sr')
	
result=import_fits(inputfile)
