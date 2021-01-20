# close all casa environments
# update pip in casa: !python -m pip install --upgrade pip
# install pandas: !python -m pip install pandas==0.22.0
# install astropy: !python -m pip install astropy

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.font_manager import FontProperties
import pandas as pd
import numpy as np
import glob
import os
import astropy.units as u
from astropy.constants import c, k_B
from astropy.wcs import WCS
from astropy.io import fits
from radio_beam import Beam
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize)




images = [
'fits/I2330P60.fits',
'smoothed_fits/45arcsec/4-5_45smoothed.fits',
'smoothed_fits/45arcsec/5-6_45smoothed.fits',
'smoothed_fits/45arcsec/6-7_45smoothed.fits',
'smoothed_fits/45arcsec/7-8_45smoothed.fits',
'smoothed_fits/45arcsec/8-9_45smoothed.fits',
'smoothed_fits/45arcsec/9-10_45smoothed.fits',
'smoothed_fits/45arcsec/10-11_45smoothed.fits',
'fits/I2330P60.fits',
]



rms = [
0.0001, # 4-5
0.00015, # 5-6
0.00013, # 6-7
0.0001, # 7-8
0.0001, # 8-9
0.0001, # 9-10
0.0001, # 10-11
0.00065, #NVSS
]



def jyb_to_mjsr(bmin,bmaj):
  fwhm_to_sigma = 1./4*np.log(2)
  beam_area = np.pi*(bmaj*bmin*fwhm_to_sigma)
  return (((1*u.Jy)/beam_area).to((1*u.MJy)/u.sr)).value


def s_nu(nu,s0,alpha): # f(xdata,a0,a1)
  return s0*nu**alpha




def em_mes(S, nu):
	T=8000 #K
	S_erg = S * 1e-17	 #this is the surface brightness in units of ergs/cm^2/sr
	c=3e10 #speed of light in cgs
	k_b=1.38e-16 #boltzmann constant in cgs
	emission_measure = (-1 *np.log(1 - ((S_erg*c**2)/(2*k_b*T*((nu*1e9)**2)))) * 1/(3.28e-7) * (T/1e4)**1.35 * (nu)**2.1)/(1e3)
	return emission_measure

def opt_dep(nu, S):
	m=em_mes(S, nu)
	t0=8e3
	tau=3.28e-7 * (t0/1e4)**(-1.35) * (nu)**(-2.1) * m
	return tau


region = open('tau_region.reg', 'r')


fluxes = []
freqs  = []
flux_err   = []
em     = []
tau	 = []

for im in images:

  header = fits.getheader(im)


  if header['OBSERVER'] == 'NVSS GRP':
    bmaj = 45.*u.arcsec
    bmin = 45.*u.arcsec
  else:
    bmaj  = imhead(im,mode='get',hdkey='beammajor')['value']*u.arcsec
    bmin  = imhead(im,mode='get',hdkey='beamminor')['value']*u.arcsec
  equiv = jyb_to_mjsr(bmin,bmaj)



  fluxes.append(imstat(im,region=region)['mean'][0]*equiv) #box='1250,1290,1290,1310'
  flux_err.append(rms[images.index(im)]*equiv)
  if header['CTYPE3'] == 'FREQ':
    freqs.append(header['CRVAL3']/1e9)
  elif header['CTYPE4'] == 'FREQ':
    freqs.append(header['CRVAL4']/1e9)
  em.append(em_mes(freqs[-1],fluxes[-1]))
  tau.append(opt_dep(freqs[-1], fluxes[-1]))



data_all = pd.DataFrame({'Frequencies (GHz)':np.round(freqs,3), 'Fluxes (MJy sr-1)':np.round(fluxes,3),
                           'RMS (Jy sr-1)': np.round(yerr,3) , 'EM (pc cm-6)': np.round(em,3),
						   'Optical Depth': np.round(tau, 3)},
                           columns=['Frequencies (GHz)', 'Fluxes (MJy sr-1)', 'RMS (Jy sr-1)', 'EM (pc cm-6)'])




ar=np.array(flux_err)
ar2=np.array(freqs)
yerr=opt_dep(ar2, ar)
xerr=0.5



fig = plt.figure()
ax = fig.add_subplot(111)
#ax.set_xscale("log");
ax.set_yscale("log")
ax.scatter(freqs, tau)

ax.errorbar(freqs,tau,
yerr=yerr,
xerr=0.5,
fmt='r+',label='data')

ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Optical Depth')
#plt.legend(loc='lower left')
plt.show()
