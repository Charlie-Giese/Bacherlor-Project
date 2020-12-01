

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



def args(argv):
	arguments=['','','','']
	try:
		opts, args = getopt.getopt(argv,"hi:p:c:o:")
	except 	getopt.GetoptError:
		print('casa_plot.py -i <inputfile> -p <power_scale_factor> -c <contour_bias_factor> -o <outputfile')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('Usage:')
			print('-i <inputfile>, required. Path to fits file relative to current working directory.')
			print('-p <scale_power_factor>, required. 0 for linear scaling, negative value to make make faint features brighter and vice versa')
			print('-c <contour_bias_factor>, optional. If blank, no contours. Takes values in range [0-inf]. Value of 0 results in contours biased towards maximum of data range, higher values scale towards min')
			print('-o <outputfile>, optional. Name of output png file, if blank does not save image')
			sys.exit()
		elif opt in ('-i'):
			inputfile = arg
			arguments[0] = arg
		elif opt in ('-p'):
			p = arg
			arguments[1] = arg
		elif opt in ('-c'):
			c = arg
			arguments[2] = arg
		elif opt in ('-o'):
			outputfile = arg
			arguments[3] = outputfile
	return arguments


def import_fits(filename):
	"""Imports a FITS radio image from CASA pipeline """
	
	fits_image_filename = os.getcwd()+filename
	with fits.open(fits_image_filename) as hdul:
		#print(hdul.info())
		data=hdul[0].data

		
	
	#now convert to Jy/str surface brigthness
	beam_fwhm = 7.2*u.arcsec
	fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
	beam_sigma = beam_fwhm * fwhm_to_sigma
	omega_B = 2 * np.pi * beam_sigma**2
	surface_brightness=((data*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(omega_B)))[0,0,:,:]
	print(np.max(surface_brightness))
	return surface_brightness


def power_scale(b_map, p, d_range):
	"""Takes as inputs a surface brightness map in Jy/str, a value for the power scaling and a data_range (list)
	   First, data is clipped so it only lies within the given data range. It is then scaled to lie between 0 and
	   and 10^p where p is the scale. Then take log of these values, giving a range of 1-p. These values will be 
	   scaled linearly to the set of availabel colors. (NOTE: This assumes negative values of p, can also take 
	   positive values, not described here"""
	
	data=b_map.value #[0,0,:,:] only keep dimensions we want


	if p < 0:
		power=-1.*p
		scaled_data=((data-np.min(data))/np.max(data))*(10**power)
		#print(np.min(scaled_data), np.max(scaled_data), scaled_data)
		output=np.log10(scaled_data+1)
		#print(np.min(output), np.max(output), output)
		return output
	elif p >0:
		scaled_data=((data-np.min(data))/np.max(data))*(p)
		output=10**scaled_data
		return output
	else:
		return data

		
def image_plot(image_data, contours, c, savefile, outputfile):
	"""Takes a scaled image in the form of a numpy array and plots it along with coordinates and a colorbar
	   The contours option takes a boolean"""
	
	
	hdu = fits.open(os.getcwd()+inputfile)
	hrd=hdu[0].header
	wcs = WCS(hrd)
	
	data=image_data.value
	
	norm = ImageNormalize(data, interval=MinMaxInterval(),
                      stretch=SqrtStretch())
	
	figure=plt.figure(num=1)
	ax=figure.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	main_image=ax.imshow(X=data, cmap='plasma', origin='lower', norm=norm, vmax=np.max(data)- 5 , vmin=np.min(data))
	cbar=figure.colorbar(main_image)
	ax.set_xlabel('Right Ascension J2000')
	ax.set_ylabel('Declination J2000')
	cbar.set_label('Surface Brigthness (MJy/Sr)')
	if contours==True:
		n = 7 #number of contours
		ub = np.max(data) #- (p-1)/(n)
		lb = np.min(data) #+ (p-1)/(n)
		spacing = c
		def level_func(lb, ub, n, spacing=1.1):
			span = (ub-lb)
			dx = 1.0 / (n-1)
			return [lb + (i*dx)**spacing*span for i in range(n)]
		levels=level_func(lb, ub, n, spacing=float(c))
		ax.contour(data, levels=levels, colors='white', alpha=0.5)
	plt.show()
	if savefile == True:
		plt.savefig(os.getcwd()+outputfile)
		

def emission_measure(S):
	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = S * 10**-23 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = -1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1	
	
	hdu = fits.open(os.getcwd()+inputfile)
	hrd=hdu[0].header
	wcs = WCS(hrd)
	
	fig=plt.figure(1) 
	ax=fig.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	em_map=ax.imshow(emission_measure, origin='lower')
	ax.set_xlabel('Right Ascension')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Emission Measure')
	plt.show()
	
	return emission_measure


arguments = args(sys.argv[1:])

inputfile = arguments[0]
p = int(arguments[1])
c = arguments[2]
outputfile = arguments[3]

if c != '':
	contours=True
else:
	contours=False

if outputfile != '':
	savefile=True
else:
	savefile=False

b_map=import_fits(inputfile)
lower_lim=np.min(b_map)
upper_lim=np.max(b_map) 
#image=power_scale(b_map, p, [lower_lim, upper_lim]) 
c=image_plot(b_map, contours, c, savefile, outputfile)
em=emission_measure(b_map.value)








