import numpy as np
import matplotlib.pyplot as plt
import math as m
import matplotlib
from scipy.interpolate import interp1d

matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['figure.figsize'] = [6.8,5.5]
matplotlib.rcParams['figure.dpi'] = 120
matplotlib.rcParams['font.sans-serif'] = "Nimbus Roman"



def flux(nu, alpha):
	flux = nu**alpha
	return flux

n = 100	

freq_thick = np.linspace(0.1, 0.5	, n)
freq_thin = np.linspace(1.5, 15, n)

optical_thick = flux(freq_thick, 2)
optical_thin = flux(freq_thin, -0.1)

x_interp = np.empty(shape=(2*n,))
y_interp = np.empty(shape=(2*n,))

x_interp[0:n] = freq_thick
x_interp[n:2*n] = freq_thin

y_interp[0:n] = optical_thick
y_interp[n:2*n] = optical_thin

f = interp1d(x_interp, y_interp, kind='quadratic')
x_new = np.linspace(0.1, 15, 1000)

spec_plot = plt.figure(num=1)
ax = spec_plot.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')

ax.plot(x_new, f(x_new), c='k')
ax.set_xlabel('\u03BD')
ax.set_ylabel('Flux Density')
ax.scatter(freq_thin, optical_thin, s=0.5)
ax.scatter(freq_thick, optical_thick, s=0.5)




plt.show()
