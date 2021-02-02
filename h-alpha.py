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
                                   ImageNormalize, LogStretch)
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


flux=data*equiv*42545250225.0
wl=header['PHOTPLAM'] #Angstrom
wl_um=wl * 1e-4


N_H = 3.944e21 #atoms cm−2

A_V = N_H/(1.9e21)
print(A_V)
R_V=3.1
X=1/wl_um

A_lambda = A_V * X / R_V

flux_true = flux * 10**(0.318 * A_V)

EM_ha = flux_true / 1.17e-7

wcs = WCS(header)

plot_data=EM_ha
norm = ImageNormalize(plot_data, interval=MinMaxInterval(), stretch=LogStretch())

figure=plt.figure(num=1)
ax=figure.add_subplot(111, projection=wcs, slices=('x','y'))
main_image=ax.imshow(X=plot_data, cmap='plasma', origin='lower', norm=norm, vmax=3e6 , vmin=0.)
cbar=figure.colorbar(main_image)
plt.show()
