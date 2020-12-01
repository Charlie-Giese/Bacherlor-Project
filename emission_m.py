
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
import sys, getopt
from astropy.wcs import WCS
from matplotlib import cm	


filename='/home/charlie/college/sophister_proj/fits/ngc7635_g0.5_briggs1.0_nsig5.image.tt0.fits'

def import_fits(filename):
	"""Imports a FITS radio image from CASA pipeline and converts it's units to Jy/sr """
	
	fits_image_filename = filename
	with fits.open(fits_image_filename) as hdul:
		#print(hdul.info())
		data=hdul[0].data

		#now convert to Jy/str surface brigthness
	beam_fwhm = 7.2*u.arcsec
	fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
	beam_sigma = beam_fwhm * fwhm_to_sigma
	omega_B = 2 * np.pi * beam_sigma**2
	surface_brightness=((data*u.Jy/u.beam).to(u.Jy/u.sr, equivalencies=u.beam_angular_area(omega_B)))[0,0,:,:] * (u.sr/u.Jy)
	return surface_brightness
	
def emission_measure(S):
	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = S * 10**-23 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = -1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1	
	return emission_measure




sb=import_fits(filename)
em=emission_measure(T, v, sb)

fig=plt.figure(1)
ax=fig.add_subplot(111)
em_map=ax.imshow(em, origin='lower')
fig.colorbar(em_map)
plt.show()
	
	
	
	

	
	
	
	
	
	
	
	
