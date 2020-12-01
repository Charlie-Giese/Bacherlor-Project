

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
import sys, getopt
from astropy.wcs import WCS
from matplotlib import cm	


filename='/home/charlie/college/project/fits/ngc7635_g0.5_briggs1.0_nsig5.image.tt0.fits'

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
	surface_brightness=((data*u.Jy/u.beam).to(u.Jy/u.sr, equivalencies=u.beam_angular_area(omega_B)))[0,0,:,:]
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
	
	
	hdu = fits.open(filename)
	hrd=hdu[0].header
	wcs = WCS(hrd)
	
	figure=plt.figure(num=1)
	ax=figure.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	main_image=ax.imshow(X=image_data, cmap='plasma', origin='lower')
	figure.colorbar(main_image)
	ax.set_xlabel('Right Ascension J2000')
	ax.set_ylabel('Declination J2000')
	if contours==True:
		n = 7 #number of contours
		ub = np.max(image_data) #- (p-1)/(n)
		lb = np.min(image_data) #+ (p-1)/(n)
		spacing = c
		def level_func(lb, ub, n, spacing=1.1):
			span = (ub-lb)
			dx = 1.0 / (n-1)
			return [lb + (i*dx)**spacing*span for i in range(n)]
		levels=level_func(lb, ub, n, spacing=1.1)
		norm = cm.colors.Normalize(vmax=abs(image_data).max(), vmin=-abs(image_data).max())
		ax.contour(image_data, levels=levels, colors='white', alpha=0.5)
	plt.show()
	if savefile == True:
		plt.savefig(os.getcwd()+outputfile)
	 

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
image=power_scale(b_map, p, [lower_lim, upper_lim]) 
c=image_plot(image, contours, c, savefile, outputfile)









