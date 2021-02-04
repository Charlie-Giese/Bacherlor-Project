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
from scipy.ndimage import gaussian_filter1d



fontsize=11
font = {'family' : 'DejaVu Sans',
'size' : fontsize}


# IMPORTING RADIO DATA AND MAKING NECESSARY CONVERSIONS

radio_image_filename = 'fits/ngc7635_g0.5_briggs1.0_nsig5.image.tt0.fits'

with fits.open(radio_image_filename) as hdul_1:
	data_1=hdul_1[0].data[0,0,:,:]
header_R = fits.getheader(radio_image_filename)

beam = Beam.from_fits_header(header_R)
SB_1=((data_1*u.Jy/u.beam).to(u.MJy/u.sr, equivalencies=u.beam_angular_area(beam))).value
SB_1masked=np.ma.masked_invalid(SB_1)

def em(SB_1masked):

	T=8000 #K
	v=8e9 # central frequency of our broadband observations
	S_erg = SB_1masked * 10**-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = (-1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*(v**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (v/1e9)**2.1)
	return emission_measure

em=em(SB_1masked)

wcs_R = WCS(header_R)[0,0,:,:]


#IMPORTING H-ALPHA DATA AND MAKING NECESSARY CONVERSIONS

inputfile_0='fits/hlsp_heritage_hst_wfc3-uvis_bubble_nebula_f656n_v1_drc.fits'
with fits.open(inputfile_0) as hdul:
	data_0=hdul[0].data
header_H = fits.getheader(inputfile_0)

wcs_H = WCS(header_H)

where_are_NaNs = np.isnan(data_0)
data_0[where_are_NaNs] = 0.0

equiv_0=header_H['PHOTFLAM'] * header_H['PHOTBW'] / header_H['D024ISCL']**2

flux_0=data_0*equiv_0*42545250225.0
wl_0=header_H['PHOTPLAM'] #Angstrom
wl_0_um=wl_0 * 1e-4

N_H = 3.944e21 #atoms cm−2
A_V = N_H/(1.9e21)
R_V=3.1
X=1/wl_0_um
A_lambda = A_V * X / R_V
flux_true = flux_0 * 10**(0.318 * A_V)

EM_ha = flux_true / 1.17e-7

"""SETTING UP THE COORDINATES OF THE 4 LINES"""

l_coord_1 = SkyCoord(350.22, 61.20, unit='deg', frame='fk5')
r_coord_1 = SkyCoord(350.15, 61.20, unit='deg', frame='fk5')

l_coord_2 = SkyCoord(350.22, 61.19, unit='deg', frame='fk5')
r_coord_2 = SkyCoord(350.15, 61.19, unit='deg', frame='fk5')

l_coord_3 = SkyCoord(350.22, 61.18, unit='deg', frame='fk5')
r_coord_3 = SkyCoord(350.15, 61.18, unit='deg', frame='fk5')

l_coord_4 = SkyCoord(350.22, 61.17, unit='deg', frame='fk5')
r_coord_4 = SkyCoord(350.15, 61.17, unit='deg', frame='fk5')

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

l_pixel_radio_4 = wcs_R.world_to_array_index(l_coord_4)
r_pixel_radio_4 = wcs_R.world_to_array_index(r_coord_4)
l_pixel_h_alpha_4 = wcs_H.world_to_array_index(l_coord_4)
r_pixel_h_alpha_4 = wcs_H.world_to_array_index(r_coord_4)

"""EXTRACTING EMISSION MEASURE VALUES"""

radio_em_vals_1 = em[l_pixel_radio_1[0] , l_pixel_radio_1[1] : r_pixel_radio_1[1]]
radio_em_vals_2 = em[l_pixel_radio_2[0] , l_pixel_radio_2[1] : r_pixel_radio_2[1]]
radio_em_vals_3 = em[l_pixel_radio_3[0] , l_pixel_radio_3[1] : r_pixel_radio_3[1]]
radio_em_vals_4 = em[l_pixel_radio_4[0] , l_pixel_radio_4[1] : r_pixel_radio_4[1]]

"""So along the line, there are 3064/173 more values in Halpha than radio"""
#halpha_pixel_size = header_H['D024SCAL'] #pixel size (arcsec) of output image

h_alpha_em_vals_1 = EM_ha[l_pixel_h_alpha_1[0] , l_pixel_h_alpha_1[1] : r_pixel_h_alpha_1[1]]
h_alpha_em_vals_2 = EM_ha[l_pixel_h_alpha_2[0] , l_pixel_h_alpha_2[1] : r_pixel_h_alpha_2[1]]
h_alpha_em_vals_3 = EM_ha[l_pixel_h_alpha_3[0] , l_pixel_h_alpha_3[1] : r_pixel_h_alpha_3[1]]
h_alpha_em_vals_4 = EM_ha[l_pixel_h_alpha_4[0] , l_pixel_h_alpha_4[1] : r_pixel_h_alpha_4[1]]

print(np.shape(radio_em_vals_1), np.shape(h_alpha_em_vals_1))

"""PLOTTING EMISSION MEASURE VALUES"""

x_R1 = np.linspace(0, 252, len(radio_em_vals_1))
x_H1 = np.linspace(0, 252, len(h_alpha_em_vals_1))
x_R2 = np.linspace(0, 252, len(radio_em_vals_2))
x_H2 = np.linspace(0, 252, len(h_alpha_em_vals_2))
x_R3 = np.linspace(0, 252, len(radio_em_vals_3))
x_H3 = np.linspace(0, 252, len(h_alpha_em_vals_3))
x_R4 = np.linspace(0, 252, len(radio_em_vals_4))
x_H4 = np.linspace(0, 252, len(h_alpha_em_vals_4))

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(x_R1, radio_em_vals_1)
axs[0, 0].plot(x_H1, h_alpha_em_vals_1)
axs[0, 0].set_title('Dec = 61.20')
axs[0, 1].plot(x_R2, radio_em_vals_2)
axs[0, 1].plot(x_H2, h_alpha_em_vals_2)
axs[0, 1].set_title('Dec = 61.19')
axs[1, 0].plot(x_R3, radio_em_vals_3)
axs[1, 0].plot(x_H3, h_alpha_em_vals_3)
axs[1, 0].set_title('Dec = 61.18')
axs[1, 1].plot(x_R4, radio_em_vals_4)
axs[1, 1].plot(x_H4, h_alpha_em_vals_4)
axs[1, 1].set_title('Dec = 61.17')

for ax in axs.flat:
    ax.set(xlabel='Arcseconds West of 350.22', ylabel='Emission Measure, pc cm^-6')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

plt.show()



"""PLOTTING RADIO DATA"""

norm = ImageNormalize(em, interval=MinMaxInterval(), stretch=LogStretch())
fig_R=plt.figure(2)
ax=fig_R.add_subplot(111, projection=wcs_R, slices=('x','y'))
em_map=ax.imshow(em, origin='lower', cmap='plasma', norm=norm, vmax=3e6, vmin=0)
ax.set_xlabel('Right Ascension\nJ2000')
ax.set_ylabel('Declination')
cbar=fig_R.colorbar(em_map)
cbar.set_label('Emission Measure, pc cm^-6')

ax.plot((l_pixel_radio_1[1], r_pixel_radio_1[1]), (l_pixel_radio_1[0], r_pixel_radio_1[0]), c='white')
ax.plot((l_pixel_radio_2[1], r_pixel_radio_2[1]), (l_pixel_radio_2[0], r_pixel_radio_2[0]), c='white')
ax.plot((l_pixel_radio_3[1], r_pixel_radio_3[1]), (l_pixel_radio_3[0], r_pixel_radio_3[0]), c='white')
ax.plot((l_pixel_radio_4[1], r_pixel_radio_4[1]), (l_pixel_radio_4[0], r_pixel_radio_4[0]), c='white')

dims=np.shape(em)
centre=(dims[0]/2., dims[1]/2.)
ax.set_xlim(centre[0]-300, centre[0]+300)
ax.set_ylim(centre[1]-300, centre[1]+300)

ra = ax.coords[0]
ra.set_format_unit('degree', decimal=True)

dec=ax.coords[1]
dec.set_format_unit('degree', decimal=True)
plt.show()

"""PLOTTING H-ALPHA DATA"""

plot_data=EM_ha
norm = ImageNormalize(plot_data, interval=MinMaxInterval(), stretch=LogStretch())
fig_H=plt.figure(num=3)
ax_0=fig_H.add_subplot(111, projection=wcs_H, slices=('x','y'))
main_image=ax_0.imshow(X=plot_data, cmap='plasma', origin='lower', norm=norm, vmax=3e6 , vmin=0.)
cbar=fig_H.colorbar(main_image)
ax_0.plot((l_pixel_h_alpha_1[1], r_pixel_h_alpha_1[1]), (l_pixel_h_alpha_1[0], r_pixel_h_alpha_1[0]), c='white')
ax_0.plot((l_pixel_h_alpha_2[1], r_pixel_h_alpha_2[1]), (l_pixel_h_alpha_2[0], r_pixel_h_alpha_2[0]), c='white')
ax_0.plot((l_pixel_h_alpha_3[1], r_pixel_h_alpha_3[1]), (l_pixel_h_alpha_3[0], r_pixel_h_alpha_3[0]), c='white')
ax_0.plot((l_pixel_h_alpha_4[1], r_pixel_h_alpha_4[1]), (l_pixel_h_alpha_4[0], r_pixel_h_alpha_4[0]), c='white')

ra = ax_0.coords[0]
ra.set_format_unit('degree', decimal=True)

dec=ax_0.coords[1]
dec.set_format_unit('degree', decimal=True)

plt.show()


"""PLOTTING RATIO OF BOTH"""
