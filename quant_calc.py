

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
import sys, getopt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from matplotlib import cm
from astropy.visualization import (MinMaxInterval, LogStretch,
                                   ImageNormalize, SqrtStretch)
import argparse
from radio_beam import Beam
from astropy.visualization.wcsaxes import SphericalCircle, Quadrangle
import matplotlib


parser = argparse.ArgumentParser(description='Calculate Emission measure and Optical depth as function of position')
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-i", "--inputfile", help="Path of input fits file relative to current working directory")

args = parser.parse_args()

inputfile=args.inputfile

matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['figure.figsize'] = [6.8,5.5]
matplotlib.rcParams['figure.dpi'] = 120
matplotlib.rcParams['font.sans-serif'] = "Nimbus Roman"



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
	em_map=ax.imshow(tau, origin='lower', cmap='viridis', norm=norm)
	ax.set_xlabel('Right Ascension\nJ2000')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Optical Depth')

	dims=np.shape(tau)
	centre=(dims[0]/2., dims[1]/2.)
	ax.set_xlim(centre[0]-300, centre[0]+300)
	ax.set_ylim(centre[1]-300, centre[1]+300)

	ra = ax.coords[0]
	ra.set_format_unit('degree', decimal=True)
	#ra.set_ticks(number=4)

	dec=ax.coords[1]
	dec.set_format_unit('degree', decimal=True)
	#dec.set_ticks(number=4)
	#ax.set_title('8-12 GHz, 5\u03C3 mask ')

	#ax.grid(True)

	plt.show()


def em(inputfile):
	S=import_fits(inputfile)[0]

	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = S * 10**-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = (-1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1)



	hrd = import_fits(inputfile)[1]
	wcs = WCS(hrd, naxis=2)


	norm = ImageNormalize(emission_measure, interval=MinMaxInterval(),
                      stretch=SqrtStretch())

	fig=plt.figure(1)
	fig.subplots_adjust(
		top=0.952,
		bottom=0.113,
		left=0.201,
		right=0.862,
		hspace=0.125,
		wspace=0.2
	)
	ax=fig.add_subplot(111, projection=wcs, slices=('x','y'))
	em_map=ax.imshow(emission_measure, origin='lower', cmap='viridis', norm=norm, vmax=np.max(emission_measure), vmin=np.min(emission_measure))
	ax.set_xlabel('Right Ascension\nJ2000')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Emission Measure, pc cm^-6')

	dims=np.shape(emission_measure)
	centre=(dims[0]/2., dims[1]/2.)
	ax.set_xlim(centre[0]-300, centre[0]+300)
	ax.set_ylim(centre[1]-300, centre[1]+300)

	ra = ax.coords[0]
	ra.set_format_unit('degree', decimal=True)

	dec=ax.coords[1]
	dec.set_format_unit('degree', decimal=True)

	#ax.grid(True)
	#350.19 61.22
	"""This following code is only for overplotting the geometry for the density estimate"""

	knot_coord = SkyCoord(350.1925, 61.2225, unit='deg', frame='fk5')
	star_coord = SkyCoord("23h20m44.5s +61d11m40.5s", frame = 'fk5')
	star_x = wcs.world_to_pixel(star_coord)[0]
	star_y = wcs.world_to_pixel(star_coord)[1]
	star = ax.scatter(star_x, star_y, marker='*', c='lime')
	knot_x = wcs.world_to_pixel(knot_coord)[0]
	knot_y = wcs.world_to_pixel(knot_coord)[1]
	knot = ax.scatter(knot_x, knot_y, marker='4', c='lime')
	x = [star_x, knot_x]
	y = [star_y, knot_y]
	ax.plot(x,y, c='lime')

	left_cord = SkyCoord(350.21, 61.202, unit='deg', frame='fk5')
	left_x = wcs.world_to_pixel(left_cord)[0]
	left_y = wcs.world_to_pixel(left_cord)[1]
	right_cord = SkyCoord(350.164, 61.205, unit='deg', frame='fk5')
	right_x = wcs.world_to_pixel(right_cord)[0]
	right_y = wcs.world_to_pixel(right_cord)[1]
	ax.plot((left_x, right_x), (left_y,right_y), c='lime')

	dens_coord = SkyCoord(350.188, 61.2035, unit='deg', frame='fk5')
	dens_x = wcs.world_to_array_index(dens_coord)[0]
	dens_y = wcs.world_to_array_index(dens_coord)[1]
	#ax.scatter(dens_x, dens_y, s=0.1)

	L = m.sqrt((350.21 - 350.164)**2 - (61.202 - 61.205)**2) * 60 * 60 #apparent angular chord length in this plot, arcseconds
	dL = L/10.
	d = 2993.440 # distance to nebula in pc
	dd = 136.388 #error, pc
	CL = L * d * 1/206265 * m.cos(61.2035 * m.pi/180) #pc, true chord length in pc
	dCL = CL * m.sqrt( (dd/d)**2 + (dL/L)**2 )
	dens_point = wcs.world_to_pixel(dens_coord)
	em_val = emission_measure[dens_x, dens_y]
	d_em_val = 1000.

	def density(value, CL):
		x = float(em_val)/float(CL)
		density = x**0.5 # per cubic cm:)
		return density

	density = density(em, CL)
	ddensity = 0.5 * m.sqrt( (d_em_val/em_val)**2 + (dCL/CL)**2 )
	print('Angular length of chord is:', L, 'pm', dL, '"')
	print('Length of Chord in pc:', CL, 'pm', dCL, 'pc')
	print('Emission Measure is', em_val, 'pm', d_em_val, 'pc cm^-6')
	print('Density is', density, 'pm', ddensity, 'cm^-3')

	plt.show()


region=[350.138, 61.2, 0.04, 0.02]



emission_measure=em(inputfile)
tau=optical_depth(inputfile)
