import matplotlib
matplotlib.use('Agg')
import aplpy
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = 'Nimbus Roman'
plt.rcParams['font.size'] = 12
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['figure.figsize']  = [10,10]
plt.rcParams['figure.dpi']  = 150
from astropy import units as u
import numpy as np
from astropy.io import fits
from astropy import stats

f_bd43 = 'bd43_n2_b0.5_pb4.image.fits'
d_bd43 = fits.getdata(f_bd43)[0,0]
f_bd60 = 'ngc_n2_b0.5_pb4.image.fits'

#*subplot*: [ xmin, ymin, dx, dy ] from the bottom left

fig = plt.figure(figsize=(6, 10))

f1 = aplpy.FITSFigure(f_bd43, figure=fig, subplot=[0.2,0.1,0.7,0.35])
f1.set_theme('publication')
f1.tick_labels.set_yformat('dd:mm')
f1.tick_labels.set_xformat('hh:mm')
#f1.set_tick_labels_font(size='x-small')
#f1.set_axis_labels_font(size='small')
f1.show_grayscale(vmin=0,vmax=4e-3)


# get sigma
mean, median, sigma_bd43 = stats.sigma_clipped_stats(d_bd43, sigma=2, iters=5)
print('RMS ',sigma_bd43)

ax2 = fig.add_axes([0.2,0.5,0.7,0.15])
dat = d_bd43.flatten()
ax2.plot(dat)
ax2.xaxis.set_visible(False)

'''
f2 = aplpy.FITSFigure(f_bd60, figure=fig, subplot=[0.5,0.1,0.35,0.8])
f2.set_theme('publication')
f2.tick_labels.set_yformat('dd:mm')
f2.tick_labels.set_xformat('hh:mm')
#f1.set_tick_labels_font(size='x-small')
#f1.set_axis_labels_font(size='small')
f2.show_grayscale(vmin=0,vmax=4e-2)

f2.hide_yaxis_label()
#f2.hide_ytick_labels()
'''
fig.savefig('my_ha.png')
#fig.close()
