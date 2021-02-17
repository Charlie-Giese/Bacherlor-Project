

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import math as m
import os
import sys, getopt
from astropy.wcs import WCS
from matplotlib import cm
from astropy.visualization import (MinMaxInterval, LogStretch,
                                   ImageNormalize, SqrtStretch)
import argparse
from radio_beam import Beam
from astropy.visualization.wcsaxes import SphericalCircle, Quadrangle


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
	em_map=ax.imshow(tau, origin='lower', cmap='viridis', norm=norm)
	ax.set_xlabel('Right Ascension\nJ2000')
	ax.set_ylabel('Declination')
	cbar=fig.colorbar(em_map)
	cbar.set_label('Optical Depth')

	dims=np.shape(tau)
	centre=(dims[0]/2., dims[1]/2.)
	ax.set_xlim(centre[0]-300, centre[0]+300)
	ax.set_ylim(centre[1]-300, centre[1]+300)

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
	wcs = WCS(hrd)


	norm = ImageNormalize(emission_measure, interval=MinMaxInterval(),
                      stretch=LogStretch())

	fig=plt.figure(1)
	ax=fig.add_subplot(111, projection=wcs, slices=('x','y',0,0))
	em_map=ax.imshow(emission_measure, origin='lower', cmap='plasma', norm=norm, vmax=3e6, vmin=0)
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

	"""This following code is only for overplotting the geometry for the density estimate"""

	#point=SphericalCircle((350.178, 61.2)*u.degree, 0.001*u.degree, edgecolor='white', facecolor='white',
	#					  transform=ax.get_transform('fk5'))
	#ax.add_patch(point)

	anchor_x= 350.138
	anchor_y = 61.2
	chord = Quadrangle((anchor_x, anchor_y)*u.degree, (0.08)*u.degree, 0.0*u.degree, vertex_unit='degree',
                       label='labels[i]', edgecolor='k', facecolor='none', linestyle='-',
                       transform=ax.get_transform('fk5'))
	#ax.scatter(anchor_x, anchor_y, transform=ax.get_transform('fk5'))
	ax.add_patch(chord)
	ax.text(350.168, 61.201, s='Chord', c='k', transform=ax.get_transform('fk5'))
	axis = Quadrangle((350.178, 61.18)*u.degree, 0.*u.degree, 0.04*u.degree, vertex_unit='degree',
	                  label='labels[i]', edgecolor='k', facecolor='none', linestyle='-',
	                  transform=ax.get_transform('fk5'))
	ax.add_patch(axis)
	ax.text(350.175, 61.176, s='Axis', c='k', transform=ax.get_transform('fk5'))

	plt.show()


region=[350.138, 61.2, 0.04, 0.02]


L = 0.08 * 60 * 60 #apparent angular chord length in this plot, arcseconds
dL = 0.01 * 60 * 60

d = 2993.440 # distance to nebula in pc
dd = 136.388 #error, pc

CL = L * d * 1/206265 * m.cos(61.2 * m.pi/180) #pc, true chord length in pc
dCL = CL * m.sqrt( (dd/d)**2 + (dL/L)**2 )


print('The angular chord length is:', L, 'pm', dL, 'arseconds')
print('The distane to the Nebula is:', d, 'pm', dd, 'pc')
print('Length of Chord in pc:', CL, 'pm', dCL, 'pc')

em_val = 29231.7 # pc/(cm^6)
d_em_val = 1000.0

def density(value, CL):
	x = float(em_val)/float(CL)
	density = x**0.5 # per cubic cm:)
	return density

density = density(em, CL)
ddensity = 0.5 * m.sqrt( (d_em_val/em_val)**2 + (dCL/CL)**2 )
print('Density is', density, 'pm', ddensity, 'cm^-3')

emission_measure=em(inputfile)
tau=optical_depth(inputfile)
