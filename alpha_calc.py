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




images = [
'alpha_fits/BN_4_5smoothed.fits',
'alpha_fits/BN_5_6smoothed.fits',
'alpha_fits/BN_6_7smoothed.fits',
'alpha_fits/BN_7_8smoothed.fits',
'alpha_fits/BN_8_9smoothed.fits',
'alpha_fits/BN_9_10smoothed.fits',
'alpha_fits/BN_10_11smoothed.fits'
]


rms = [
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

def em_mes(nu,s0):
  t0     = 8e3       # temperature K
  a      = 2*(k_B.cgs)*3.28*1e-7/(c.cgs)**2
  tau_nu2_T_em = a*1e4*(t0/1e4)**(-0.35)*1e18*(nu)**(-0.1)
  return (s0/tau_nu2_T_em*1e-23).value


regions = open('alpha_regions.txt', 'r')
print regions


for reg in regions.readlines():

  if '#' in reg: continue
  fluxes = []
  freqs  = []
  yerr   = []
  em     = []
  for im in images:

     header = fits.getheader(im)


     bmaj  = imhead(im,mode='get',hdkey='beammajor')['value']*u.arcsec
     bmin  = imhead(im,mode='get',hdkey='beamminor')['value']*u.arcsec
     equiv = jyb_to_mjsr(bmin,bmaj)
     try:
       fluxes.append(imstat(im,region=reg)['mean'][0]*equiv) #box='1250,1290,1290,1310'
       yerr.append(rms[images.index(im)]*equiv)
       freqs.append(header['RESTFRQ']/1e9)
       #freqs.append(imhead(imagename=im,mode='get',hdkey='crval4')['value']/1e9)
       em.append(em_mes(freqs[-1],fluxes[-1]))
     except:
       continue


  #area = (reg.split(',')[-3])

  data_all = pd.DataFrame({'Frequencies (GHz)':np.round(freqs,3), 'Fluxes (MJy sr-1)':np.round(fluxes,3),
                           'RMS (Jy sr-1)': np.round(yerr,3) , 'EM (pc cm-6)': np.round(em,3)},
                           columns=['Frequencies (GHz)', 'Fluxes (MJy sr-1)', 'RMS (Jy sr-1)', 'EM (pc cm-6)'])


  #print(area)
  #print(data_all)

  #fig=plt.figure()
  #ax=fig.add_subplot(111)
  #ax.scatter(freqs, fluxes)
  #plt.show()




  coef, cov = curve_fit(s_nu,freqs,fluxes,sigma=yerr,absolute_sigma=True)    #cf(f,xdata,ydata)
  s0        = coef[0]
  freq_err  = np.sqrt(cov[0,0]) # [[sxx sxy][syx syy]]
  alpha     = coef[1]
  alpha_err = np.sqrt(cov[1,1])

  print(alpha)


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
  ax.set_xscale("log"); ax.set_yscale("log")
  #ax.scatter(freqs, fluxes)
  ax.plot(freqs, y, c='k', label=r'$\alpha$'+'='+str(round(coef[1],2))+' +/-'+str(round(alpha_err,2)))

  ax.plot(freqs,yp,'g--')
  ax.plot(freqs,ym,'g--')

  ax.fill_between(freqs,yp,ym,facecolor='gray',alpha=0.15)

  ax.set_xlim([2,12.6])
  #ax.set_ylim([0.1, 8])


  ax.errorbar(freqs,fluxes,
  yerr=yerr_fixed,
  xerr=0.5,
  fmt='r+',label='data')

  plt.legend(loc='lower right')
  ax.set_xlabel('Frequency (GHz)')
  ax.set_ylabel('Surface Brightness (MJy/sr)')

  plt.show()




"""
  y  = s_nu(freqs,s0,alpha)

  yerr_fixed = np.array(yerr)

  yerr_fixed[yerr>=y] = y[yerr>=y]*.9999 # fix the values that are way too large

  ypm = s_nu(freqs,coef[0]+freq_err,coef[1]-alpha_err)
  ypp = s_nu(freqs,coef[0]+freq_err,coef[1]+alpha_err)
  yp  = np.maximum(ypm,ypp)

  ymp = s_nu(freqs,coef[0]-freq_err,coef[1]+alpha_err)
  ymm = s_nu(freqs,coef[0]-freq_err,coef[1]-alpha_err)
  ym  = np.minimum(ymp,ymm)


  #print('freq:',lgfr,'fluxes:',lgfl,'yerr:',yerr,'fit:',y)
  fig = plt.figure()
  ax  = plt.axes()
  ax.set_xscale("log"); ax.set_yscale("log",nonposy="clip")

  ax.errorbar(freqs,fluxes,
  yerr=yerr_fixed,
  xerr=0.5,
  fmt='r+',label='data')

  ax.plot(freqs,y,label= r'$\alpha$'+'='+str(round(coef[1],2))+' +/-'+str(round(alpha_err,2)))

  ax.plot(freqs,yp,'g--')
  ax.plot(freqs,ym,'g--')
  ax.fill_between(freqs,yp,ym,facecolor='gray',alpha=0.15)

  ax.set_xlim([4.4,11.6])
  ax.set_ylim([1e-5,1e-2])
  ax.set_ylabel(r'$F$ (Jy)')
  ax.set_xlabel(r'$\nu$ (GHz)')
  #ax.set_title()
  fp = FontProperties()
  fp.set_size('x-small')
  #plt.legend(title=area,#bbox_to_anchor=(1.007,1.14),
  	#loc='lower right',prop=fp,frameon=False)

  plt.show()
"""
