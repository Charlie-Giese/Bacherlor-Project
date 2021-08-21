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


matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['figure.figsize'] = [6.8,5.5]
matplotlib.rcParams['figure.dpi'] = 120
matplotlib.rcParams['font.sans-serif'] = "Nimbus Roman"

images = [
#'fits/I2330P60.fits',
'smoothed_fits/standard/unmasked/4-5_smoothed.fits',
'smoothed_fits/standard/unmasked/5-6_smoothed.fits',
'smoothed_fits/standard/unmasked/6-7_smoothed.fits',
'smoothed_fits/standard/unmasked/7-8_smoothed.fits',
'smoothed_fits/standard/unmasked/8-9_smoothed.fits',
'smoothed_fits/standard/unmasked/9-10_smoothed.fits',
'smoothed_fits/standard/unmasked/10-11_smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/4-5_45smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/5-6_45smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/6-7_45smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/7-8_45smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/8-9_45smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/9-10_45smoothed.fits',
#'smoothed_fits/45arcsec/unmasked/10-11_45smoothed.fits',
]


rms = [
0.00065, #NVSS
0.0001, # 4-5
0.00015, # 5-6
0.00013, # 6-7
0.0001, # 7-8
0.0001, # 8-9
0.0001, # 9-10
0.0001, # 10-11
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
	emission_measure = -1. *np.log(1. - ((S_erg * c**2.)/(2.*k_b*T*((nu*1e9)**2.)))) * 1./(3.28e-7) * (T/1e4)**1.35 * (nu)**2.1
	return emission_measure

def opt_dep(nu, em):
	t0=8e3
	tau=3.28e-7 * (t0/1e4)**(-1.35) * (nu)**(-2.1) * em
	return tau


regions = open('regions.reg', 'r')

titles = [
'Region (a)',
'Region (b)',
'Region (c)',
'Region (d)'
]

i=0

for reg in regions.readlines():

  if '#' in reg: continue
  fluxes = []
  freqs  = []
  yerr   = []
  em     = []
  tau    = []
  for im in images:


     header = fits.getheader(im)

     if header['OBSERVER'] == 'NVSS GRP':
       bmaj = 45.*u.arcsec
       bmin = 45.*u.arcsec
     else:
       bmaj  = imhead(im,mode='get',hdkey='beammajor')['value']*u.arcsec
       bmin  = imhead(im,mode='get',hdkey='beamminor')['value']*u.arcsec
     equiv = jyb_to_mjsr(bmin,bmaj)

     try:
       fluxes.append(imstat(im,region=reg)['mean'][0]*equiv) #box='1250,1290,1290,1310'
       yerr.append(rms[images.index(im)]*equiv)
       if header['CTYPE3'] == 'FREQ':
         freqs.append(header['CRVAL3']/1e9)
       elif header['CTYPE4'] == 'FREQ':
         freqs.append(header['CRVAL4']/1e9)
       #freqs.append(imhead(imagename=im,mode='get',hdkey='crval4')['value']/1e9)
       em.append(em_mes(freqs[-1],fluxes[-1]))
       tau.append(opt_dep(freqs[-1], em_mes(freqs[-1],fluxes[-1])))
     except:
       continue





  #data_all = pd.DataFrame({'Frequencies (GHz)':np.round(freqs,3), 'Fluxes (MJy sr-1)':np.round(fluxes,3),
    #                       'RMS (Jy sr-1)': np.round(yerr,3) , 'EM (pc cm-6)': np.round(em,3)},
     #                      columns=['Frequencies (GHz)', 'Fluxes (MJy sr-1)', 'RMS (Jy sr-1)', 'EM (pc cm-6)'])





  coef, cov = curve_fit(s_nu,freqs,fluxes,sigma=yerr,absolute_sigma=True)    #cf(f,xdata,ydata)
  s0        = coef[0]
  freq_err  = np.sqrt(cov[0,0]) # [[sxx sxy][syx syy]]
  alpha     = coef[1]
  alpha_err = np.sqrt(cov[1,1])

  print 'Spectral Index is:', alpha


  y  = s_nu(freqs,s0,alpha)
  yerr_fixed = np.array(yerr) * 3


  ypm = s_nu(freqs,coef[0]+freq_err,coef[1]-alpha_err)
  ypp = s_nu(freqs,coef[0]+freq_err,coef[1]+alpha_err)
  yp  = np.maximum(ypm,ypp)

  ymp = s_nu(freqs,coef[0]-freq_err,coef[1]+alpha_err)
  ymm = s_nu(freqs,coef[0]-freq_err,coef[1]-alpha_err)
  ym  = np.minimum(ymp,ymm)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  #ax.set_xscale("log"); ax.set_yscale("log")
  #ax.scatter(freqs, fluxes)

  ax.plot(freqs, y, c='k', label=r'$\alpha$'+'='+str(round(coef[1],2))+' +/-'+str(round(alpha_err,3)))

  ax.plot(freqs,yp,'g--')
  ax.plot(freqs,ym,'g--')

  #ax.set_yticks( [1e-2, 1e-1, 1e0, 1e1] )

  ax.fill_between(freqs,yp,ym,facecolor='gray',alpha=0.15)

  #ax.set_xlim([3.5,10**1.1])
  #ax.set_ylim([10, 16])


  ax.errorbar(freqs,fluxes,
  yerr=yerr_fixed,
  xerr=0.5,
  fmt='r+',label='data')

  ax.tick_params(axis='x', labelsize= 14)
  ax.tick_params(axis='y', labelsize= 14)

  plt.legend(loc='upper right')
  ax.set_xlabel('Frequency (GHz)')
  ax.set_ylabel('Surface Brightness (MJy/sr)')

  ax.set_title(titles[i])
  i+=1

  plt.show()
