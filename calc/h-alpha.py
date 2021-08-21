#Author: Charlie Giese
#Code for comparing VLA radio data with H-alpha image of Bubble Nebula

#Imports

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
from astropy.coordinates import SkyCoord
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
import matplotlib

#setting graphical parameters

matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['figure.figsize'] = [6.8,5.5]
matplotlib.rcParams['figure.dpi'] = 120
matplotlib.rcParams['font.sans-serif'] = "Nimbus Roman"


# IMPORTING RADIO DATA AND MAKING NECESSARY CONVERSIONS

def import_radio(file):

	print('Importing Radio band data')

	radio_image_filename = file

	with fits.open(radio_image_filename) as hdul_1:
		data_1=hdul_1[0].data[0,0,:,:]
		header_R = fits.getheader(radio_image_filename)

	beam = Beam.from_fits_header(header_R)
	SB_1masked=np.ma.masked_invalid(data_1)
	return SB_1masked


def em(radio_flux):

	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = SB_1masked * 10**-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = (-1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1)
	return emission_measure

file = sys.argv[1]
radio_flux = import_radio(file)
em=em(radio_flux)
wcs_R = WCS(header_R)[0,0,:,:]


#IMPORTING H-ALPHA DATA AND MAKING NECESSARY CONVERSIONS

print('Importing H-Alpha data')
inputfile_H= 'hlsp_heritage_hst_wfc3-uvis_bubble_nebula_f656n_v1_drc.fits'
with fits.open(inputfile_H) as hdul:
	data_H=hdul[0].data
header_H = fits.getheader(inputfile_H)

wcs_H = WCS(header_H)

where_are_NaNs = np.isnan(data_H)
data_H[where_are_NaNs] = 0.0

equiv_H=header_H['PHOTFLAM'] * header_H['PHOTBW'] / header_H['D024ISCL']**2

flux_H=data_H*equiv_H*42545250225.0
wl_H=header_H['PHOTPLAM'] #Angstrom
wl_H_um=wl_H * 1e-4

N_H=6.41e21  #atoms cm−2
A_V = N_H/(1.9e21)
R_V=3.1
X=1/wl_0_um
A_lambda = A_V * X / R_V
flux_true = flux_H * 10**(0.318 * A_V)

bmaj = 20
bmin = 20

sigma_1 = bmaj /(2 * m.sqrt(2 * m.log(2)) * header_H['D024ISCL'])
sigma_2 = bmin /(2 * m.sqrt(2 * m.log(2)) * header_H['D024ISCL'])

smoothed = gaussian_filter(flux_true, sigma = [sigma_1,sigma_2])

EM_ha = smoothed / 1.17e-7

"""SETTING UP THE COORDINATES OF THE 4 LINES"""

print('Calculating coordinates of slices')

l_coord_1 = SkyCoord(350.22, 61.202, unit='deg', frame='fk5')
r_coord_1 = SkyCoord(350.15, 61.202, unit='deg', frame='fk5')

l_coord_2 = SkyCoord(350.22, 61.197, unit='deg', frame='fk5')
r_coord_2 = SkyCoord(350.15, 61.197, unit='deg', frame='fk5')

l_coord_3 = SkyCoord(350.22, 61.192, unit='deg', frame='fk5')
r_coord_3 = SkyCoord(350.15, 61.192, unit='deg', frame='fk5')

l_pixel_radio_1 = wcs_R.world_to_array_index(l_coord_1)
r_pixel_radio_1 = wcs_R.world_to_array_index(r_coord_1)
l_pixel_h_alpha_1 = wcs_H.world_to_array_index(l_coord_1)
r_pixel_h_alpha_1 = wcs_H.world_to_array_index(r_coord_1)

l_pixel_radio_2 = wcs_R.world_to_array_index(l_coord_2)
r_pixel_radio_2 = wcs_R.world_to_array_index(r_coord_2)
l_pixel_h_alpha_2 = wcs_H.world_to_array_index(l_coord_2)
r_pixel_h_alpha_2 = wcs_H.world_to_array_index(r_coord_2)

l_pixel_radio_3 = wcs_R.world_to_array_index(l_coord_3)
r_pixel_radio_3 = wcs_R.world_to_array_index(r_coord_3)
l_pixel_h_alpha_3 = wcs_H.world_to_array_index(l_coord_3)
r_pixel_h_alpha_3 = wcs_H.world_to_array_index(r_coord_3)

"""EXTRACTING EMISSION MEASURE VALUES"""

radio_em_vals_1 = em[l_pixel_radio_1[0] , l_pixel_radio_1[1] : r_pixel_radio_1[1]]
radio_em_vals_2 = em[l_pixel_radio_2[0] , l_pixel_radio_2[1] : r_pixel_radio_2[1]]
radio_em_vals_3 = em[l_pixel_radio_3[0] , l_pixel_radio_3[1] : r_pixel_radio_3[1]]

"""So along the line, there are 3064/173 more values in Halpha than radio"""

h_alpha_em_vals_1 = EM_ha[l_pixel_h_alpha_1[0] , l_pixel_h_alpha_1[1] : r_pixel_h_alpha_1[1]]
h_alpha_em_vals_2 = EM_ha[l_pixel_h_alpha_2[0] , l_pixel_h_alpha_2[1] : r_pixel_h_alpha_2[1]]
h_alpha_em_vals_3 = EM_ha[l_pixel_h_alpha_3[0] , l_pixel_h_alpha_3[1] : r_pixel_h_alpha_3[1]]

"""SETTING UP FIGURE"""

fig = plt.figure(1)
ax1 = fig.add_subplot(311, label = 'em_comp')
ax2 = fig.add_subplot(312, label = 'radio_image', projection = wcs_R, slices=('x', 'y'))
ax3 = fig.add_subplot(313, label = 'halpha_iamge')

"""PLOTTING EMISSION MEASURE VALUES"""

x_R1 = np.linspace(0, 252, len(radio_em_vals_1))
x_H1 = np.linspace(0, 252, len(h_alpha_em_vals_1))
x_R2 = np.linspace(0, 252, len(radio_em_vals_2))
x_H2 = np.linspace(0, 252, len(h_alpha_em_vals_2))
x_R3 = np.linspace(0, 252, len(radio_em_vals_3))
x_H3 = np.linspace(0, 252, len(h_alpha_em_vals_3))

#ax1.spines['top'].set_color('none')
#ax1.spines['bottom'].set_color('none')
#ax1.spines['left'].set_color('none')
#ax1.spines['right'].set_color('none')
#ax1.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
# Set common labels
#ax1.set_xlabel('Arcseconds West of 350.22\N{DEGREE SIGN}')
#ax1.set_ylabel('Emission Measure, $pc\:cm^{-6}$', labelpad=25.)

#ax0.plot(x_R1, radio_em_vals_1)
#ax0.plot(x_H1, h_alpha_em_vals_1)
#ax0.set_yscale('log')
#ax0.set_title('Dec = 61.202\N{DEGREE SIGN}')
#ax1 = fig.add_subplot(312)
ax1.plot(x_R2, radio_em_vals_2)
ax1.plot(x_H2, h_alpha_em_vals_2)
ax1.set_yscale('log')
ax1.set_title('Dec = 61.197\N{DEGREE SIGN}')
#ax2 = fig.add_subplot(313)
#ax2.plot(x_R3, radio_em_vals_3)
#ax2.plot(x_H3, h_alpha_em_vals_3)
#ax2.set_yscale('log')
#ax2.set_title('Dec = 61.192\N{DEGREE SIGN}')


# Hide x labels and tick labels for top plots and y ticks for right plots.

#ax0.label_outer()
#ax1.label_outer()
#ax2.label_outer()

#plt.savefig(radio_image_filename+'__EM_COMP.png')
#plt.show()


"""PLOTTING RADIO EMISSION MEASURE"""

norm = ImageNormalize(em, interval=MinMaxInterval(), stretch=SqrtStretch())
#fig_R=plt.figure(2)
#ax=fig_R.add_subplot(111, projection=wcs_R, slices=('x','y'))
em_map=ax2.imshow(em, origin='lower', cmap='plasma', norm=norm, vmax=3e6, vmin=0)
ax2.set_xlabel('Right Ascension\nJ2000')
ax2.set_ylabel('Declination')
fig.colorbar(em_map, ax2)
cbar.set_label('Emission Measure, $pc\:cm^{-6}$')

ax2.plot((l_pixel_radio_1[1], r_pixel_radio_1[1]), (l_pixel_radio_1[0], r_pixel_radio_1[0]), c='white')
ax2.plot((l_pixel_radio_2[1], r_pixel_radio_2[1]), (l_pixel_radio_2[0], r_pixel_radio_2[0]), c='white')
ax2.plot((l_pixel_radio_3[1], r_pixel_radio_3[1]), (l_pixel_radio_3[0], r_pixel_radio_3[0]), c='white')

dims=np.shape(em)
centre=(dims[0]/2., dims[1]/2.)
ax2.set_xlim(centre[0]-300, centre[0]+300)
ax2.set_ylim(centre[1]-300, centre[1]+300)

ra = ax2.coords[0]
ra.set_format_unit('degree', decimal=True)
ax2.set_xlabel('Right Ascension J2000')
dec=ax2.coords[1]
dec.set_format_unit('degree', decimal=True)
ax2.set_ylabel('Declination')
#plt.savefig(radio_image_filename+'__R-EM.png')
#plt.show()


"""PLOTTING RATIO OF BOTH"""

x = np.linspace(0, 252, num = 170)

R_inter_1 = interp1d(x_R1, radio_em_vals_1, kind = 'cubic')
R_inter_2 = interp1d(x_R2, radio_em_vals_2, kind = 'cubic')
R_inter_3 = interp1d(x_R3, radio_em_vals_3, kind = 'cubic')

H_inter_1 = interp1d(x_H1, h_alpha_em_vals_1, kind = 'cubic')
H_inter_2 = interp1d(x_H2, h_alpha_em_vals_2, kind = 'cubic')
H_inter_3 = interp1d(x_H3, h_alpha_em_vals_3, kind = 'cubic')

R1 = R_inter_1(x)
R2 = R_inter_2(x)
R3 = R_inter_3(x)

H1 = H_inter_1(x)
H2 = H_inter_2(x)
H3 = H_inter_3(x)

ratio_1 = H1 / R1
ratio_2 = H2 / R2
ratio_3 = H3 / R3

ax3.plot(x, ratio_1, label='61.202\N{DEGREE SIGN}')
ax3.plot(x, ratio_2, label='61.197\N{DEGREE SIGN}')
ax3.plot(x, ratio_3, label='61.192\N{DEGREE SIGN}')

ax3.legend()
ax3.set_xlabel('Arcseconds left of 350.22\N{DEGREE SIGN}')
ax3.set_ylabel('H\u03B1 EM / Radio EM')
ax3.set_yscale('log')
#ax3.set_ylim(0, 3.)
#plt.savefig(radio_image_filename+'__RATIO.png')
plt.show()


"""PLOTTING H-ALPHA emission measure"""
"""
plot_data=EM_ha
norm = ImageNormalize(plot_data, interval=MinMaxInterval(), stretch=SqrtStretch())
fig_H1=plt.figure(num=3)
ax_01=fig_H1.add_subplot(111, projection=wcs_H, slices=('x','y'))
main_image=ax_01.imshow(X=plot_data, cmap='plasma', origin='lower', norm=norm, vmax=3e6 , vmin=0.)
cbar=fig_H1.colorbar(main_image)
cbar.set_label('Emission Measure, $pc\:cm^{-6}$')
ax_01.plot((l_pixel_h_alpha_1[1], r_pixel_h_alpha_1[1]), (l_pixel_h_alpha_1[0], r_pixel_h_alpha_1[0]), c='white')
ax_01.plot((l_pixel_h_alpha_2[1], r_pixel_h_alpha_2[1]), (l_pixel_h_alpha_2[0], r_pixel_h_alpha_2[0]), c='white')
#ax_01.plot((l_pixel_h_alpha_3[1], r_pixel_h_alpha_3[1]), (l_pixel_h_alpha_3[0], r_pixel_h_alpha_3[0]), c='white')
#ax_01.plot((l_pixel_h_alpha_4[1], r_pixel_h_alpha_4[1]), (l_pixel_h_alpha_4[0], r_pixel_h_alpha_4[0]), c='white')

ra = ax_01.coords[0]
ra.set_format_unit('degree', decimal=True)
ax_01.set_xlabel('Right Ascension J2000')
dec=ax_01.coords[1]
dec.set_format_unit('degree', decimal=True)
ax_01.set_ylabel('Declination')
#plt.savefig(radio_image_filename+'__HA-EM.png')
plt.show()
"""
#Plotting non-smoothed H-Alpha

"""
norm = ImageNormalize(flux_true, interval=MinMaxInterval(), stretch=LogStretch())
fig_H2=plt.figure(num=4)
ax_02=fig_H2.add_subplot(111, projection=wcs_H, slices=('x','y'))
main_image=ax_02.imshow(X=flux_true, cmap='plasma', origin='lower', norm=norm, vmax=3e6 , vmin=0.)
cbar=fig_H2.colorbar(main_image)
ax_01.plot((l_pixel_h_alpha_1[1], r_pixel_h_alpha_1[1]), (l_pixel_h_alpha_1[0], r_pixel_h_alpha_1[0]), c='white')
ax_01.plot((l_pixel_h_alpha_2[1], r_pixel_h_alpha_2[1]), (l_pixel_h_alpha_2[0], r_pixel_h_alpha_2[0]), c='white')
ax_01.plot((l_pixel_h_alpha_3[1], r_pixel_h_alpha_3[1]), (l_pixel_h_alpha_3[0], r_pixel_h_alpha_3[0]), c='white')
ax_01.plot((l_pixel_h_alpha_4[1], r_pixel_h_alpha_4[1]), (l_pixel_h_alpha_4[0], r_pixel_h_alpha_4[0]), c='white')

ra = ax_02.coords[0]
ra.set_format_unit('degree', decimal=True)
ax_02.set_xlabel('Right Ascension J2000')
dec=ax_02.coords[1]
dec.set_format_unit('degree', decimal=True)
ax_02.set_ylabel('Declination')
plt.savefig(radio_image_filename+'__NONSMOOTHED.png')

"""
