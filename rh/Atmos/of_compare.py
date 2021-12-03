import numpy as np
import matplotlib.pyplot as plt

rh = np.loadtxt('rh_of').T
sasha = np.loadtxt('of_sasha.dat', skiprows=1).T
lam = np.linspace(1500, 4000, num=1301)

def line_of(lam):
	lb_coeff = 2.1177*np.exp(-(lam-2087.7)*(lam-2087.7)/2.421E6) + 0.68738;
	return (1+2/3*lb_coeff)

plt.plot(rh[0], rh[1]+1, label=r'RH: H$^-$')
plt.plot(rh[0], rh[3]+1, label=r'RH: metals')
plt.plot(sasha[0]/10, sasha[1], label='Shapiro+10')
plt.plot(lam/10, line_of(lam), label='Busa+01')

plt.xlim([150, 400])
plt.xlabel(r'$\lambda$ [nm]')
plt.ylabel(r'$f(\lambda)$')

plt.legend(loc=1, frameon=False, ncol=2)
plt.savefig('compare_of.png')

plt.show()