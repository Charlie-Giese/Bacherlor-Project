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
import aplpy
import argparse
from radio_beam import Beam
from astropy.coordinates import SkyCoord
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
import matplotlib

#setting graphical parameters

matplotlib.rcParams['font.sans-serif'] = 'Nimbus Roman'
plt.rcParams['font.size'] = 12
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize']  = [10,10]
plt.rcParams['figure.dpi']  = 150


# IMPORTING RADIO DATA AND MAKING NECESSARY CONVERSIONS

def import_radio(file):

	print('Importing Radio band data')

	radio_image_filename = file

	with fits.open(radio_image_filename) as hdul_1:
		data_1=hdul_1[0].data[0,0,:,:]
		header_R = fits.getheader(radio_image_filename)

	beam = Beam.from_fits_header(header_R)
	SB_1masked=np.ma.masked_invalid(data_1)
	return SB_1masked, header_R


def em(radio_flux):

	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = radio_flux * 10**-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = (-1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1)
	return emission_measure

file = sys.argv[1]
radio_flux, header_R = import_radio(file)
em=em(radio_flux)
wcs_R = WCS(header_R)[0,0,:,:]


#IMPORTING H-ALPHA DATA AND MAKING NECESSARY CONVERSIONS
"""
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
X=1/wl_H_um
A_lambda = A_V * X / R_V
flux_true = flux_H * 10**(0.318 * A_V)

bmaj = 20
bmin = 20

sigma_1 = bmaj /(2 * m.sqrt(2 * m.log(2)) * header_H['D024ISCL'])
sigma_2 = bmin /(2 * m.sqrt(2 * m.log(2)) * header_H['D024ISCL'])

smoothed = gaussian_filter(flux_true, sigma = [sigma_1,sigma_2])

EM_ha = smoothed / 1.17e-7

#SETTING UP THE COORDINATES OF THE 4 LINES

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
"""
"""EXTRACTING EMISSION MEASURE VALUES"""
"""
radio_em_vals_1 = em[l_pixel_radio_1[0] , l_pixel_radio_1[1] : r_pixel_radio_1[1]]
radio_em_vals_2 = em[l_pixel_radio_2[0] , l_pixel_radio_2[1] : r_pixel_radio_2[1]]
radio_em_vals_3 = em[l_pixel_radio_3[0] , l_pixel_radio_3[1] : r_pixel_radio_3[1]]

#So along the line, there are 3064/173 more values in Halpha than radio

h_alpha_em_vals_1 = EM_ha[l_pixel_h_alpha_1[0] , l_pixel_h_alpha_1[1] : r_pixel_h_alpha_1[1]]
h_alpha_em_vals_2 = EM_ha[l_pixel_h_alpha_2[0] , l_pixel_h_alpha_2[1] : r_pixel_h_alpha_2[1]]
h_alpha_em_vals_3 = EM_ha[l_pixel_h_alpha_3[0] , l_pixel_h_alpha_3[1] : r_pixel_h_alpha_3[1]]
"""
"""SETTING UP FIGURE"""

fig = plt.figure(figsize=(6, 10))

f1 = aplpy.FITSFigure(em, figure=fig, subplot=[0.15,0.4,0.7,0.35])
ax2 = fig.add_axes([0.15, 0.1, 0.7, 0.25])
ax3 = fig.add_axes([0.15, 0.8, 0.7, 0.15])

plt.show()
# Set common labels for axsTop
#ax1 = subfigs[0].add_subplot(111)
#ax1.set_xlabel('Arcseconds West of 350.22\N{DEGREE SIGN}')
#ax1.set_ylabel('Emission Measure, $pc\:cm^{-6}$', labelpad=25.)
#ax1.spines['top'].set_color('none')
#ax1.spines['bottom'].set_color('none')
#ax1.spines['left'].set_color('none')
#ax1.spines['right'].set_color('none')
#ax1.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
"""
x_arrays = [
np.linspace(0, 252, len(radio_em_vals_1)),
np.linspace(0, 252, len(h_alpha_em_vals_1)),
np.linspace(0, 252, len(radio_em_vals_2)),
np.linspace(0, 252, len(h_alpha_em_vals_2)),
np.linspace(0, 252, len(radio_em_vals_3)),
np.linspace(0, 252, len(h_alpha_em_vals_3))
]
xval_arrays=[radio_em_vals_1,h_alpha_em_vals_1,radio_em_vals_2,h_alpha_em_vals_2,radio_em_vals_3,h_alpha_em_vals_3]
i,j = 0,1

#titles = ['Dec = 61.202\N{DEGREE SIGN}', 'Dec = 61.197\N{DEGREE SIGN}', 'Dec = 61.192\N{DEGREE SIGN}']
for ax in axsTop:
	ax.set_yscale('log')
	#ax.set_ylabel('Emission Measure, $pc\:cm^{-6}$', labelpad=25.)
	ax.set_xticks(ticks=[])
	ax.plot(x_arrays[i], xval_arrays[i])
	ax.plot(x_arrays[j], xval_arrays[j])
	#ax.set_title(titles[i])
	i += 1
	j +=1

axMid.legend()
axMid.set_ylabel('H\u03B1 EM / Radio EM')
axMid.set_yscale('log')
axMid.set_ylim(0, 3.)
x = np.linspace(0, 252, num = 170)

R_inter_1 = interp1d(x_arrays[0], xval_arrays[0], kind = 'cubic')
R_inter_2 = interp1d(x_arrays[2], xval_arrays[2], kind = 'cubic')
R_inter_3 = interp1d(x_arrays[4], xval_arrays[4], kind = 'cubic')

H_inter_1 = interp1d(x_arrays[1], xval_arrays[1], kind = 'cubic')
H_inter_2 = interp1d(x_arrays[3], xval_arrays[3], kind = 'cubic')
H_inter_3 = interp1d(x_arrays[5], xval_arrays[5], kind = 'cubic')

R1 = R_inter_1(x)
R2 = R_inter_2(x)
R3 = R_inter_3(x)
H1 = H_inter_1(x)
H2 = H_inter_2(x)
H3 = H_inter_3(x)
ratio_1 = H1 / R1
ratio_2 = H2 / R2
ratio_3 = H3 / R3

axMid.plot(x, ratio_1, label='61.202\N{DEGREE SIGN}')
axMid.plot(x, ratio_2, label='61.197\N{DEGREE SIGN}')
axMid.plot(x, ratio_3, label='61.192\N{DEGREE SIGN}')

axBot.set_xlabel('Right Ascension')
axBot.set_ylabel('Declination')
norm = ImageNormalize(em, interval=MinMaxInterval(), stretch=SqrtStretch())
em_map=axBot.imshow(em, origin='lower', cmap='binary', norm=norm, vmax=np.max(em), vmin=np.min(em))
cbar = fig.colorbar(em_map, axBot)
cbar.set_label('Emission Measure, $pc\:cm^{-6}$')
axBot.plot((l_pixel_radio_1[1], r_pixel_radio_1[1]), (l_pixel_radio_1[0], r_pixel_radio_1[0]), c='white')
axBot.plot((l_pixel_radio_2[1], r_pixel_radio_2[1]), (l_pixel_radio_2[0], r_pixel_radio_2[0]), c='white')
axBot.plot((l_pixel_radio_3[1], r_pixel_radio_3[1]), (l_pixel_radio_3[0], r_pixel_radio_3[0]), c='white')
dims=np.shape(em)
centre=(dims[0]/2., dims[1]/2.)
axBot.set_xlim(centre[0]-150, centre[0]+150)
axBot.set_ylim(centre[1]-150, centre[1]+150)
#ra = axBot.coords[0]
#ra.set_format_unit('degree', decimal=True)
#dec=axBot.coords[1]
#dec.set_format_unit('degree', decimal=True)
plt.show()

"""
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
