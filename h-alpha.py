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

fontsize=11
font = {'family' : 'DejaVu Sans',
'size' : fontsize}

with fits.open(inputfile) as hdul:
	data=hdul[0].data
header = fits.getheader(inputfile)

where_are_NaNs = np.isnan(data)
data[where_are_NaNs] = 0.0

equiv=header['PHOTFLAM'] * header['PHOTBW'] / header['D024ISCL']**2

print(equiv)

flux=data*equiv
wl=header['PHOTPLAM']

N_H = 2.7 × 10**21 #atoms cm−2

A_V = N_H/(1.9e21)
R_V=3.1
X=1e-6/wl

A_lambda = A_V * R_V / X



wcs = WCS(header)

norm = ImageNormalize(flux, interval=MinMaxInterval(), stretch=SqrtStretch())

figure=plt.figure(num=1)
ax=figure.add_subplot(111, projection=wcs, slices=('x','y'))
main_image=ax.imshow(X=flux, cmap='plasma', origin='lower', norm=norm, vmax=np.max(flux) , vmin=np.min(flux))
cbar=figure.colorbar(main_image)
plt.show()




def em_mes(S, nu):
	T=8000 #K
	S_erg = S * 1e-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = -1. *np.log(1. - ((S_erg * c**2.)/(2.*k_b*T*((nu*1e9)**2.)))) * 1./(3.28e-7) * (T/1e4)**1.35 * (nu)**2.1
	return emission_measure
