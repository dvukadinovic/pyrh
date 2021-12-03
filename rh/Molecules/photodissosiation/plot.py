import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from scipy.constants import h

lam1, pd1 = np.loadtxt("CH.txt", unpack=True, usecols=(0,2))
lam3, pd3 = np.loadtxt("ch_new.pdf", unpack=True, usecols=(1,2))

# plt.plot(lam1, np.log10(pd1))
# plt.plot(lam3/10, np.log10(pd3))
# plt.show()

# data.shape = 105 x 15
data = np.loadtxt("CH_kurucz.txt", delimiter=",")

temp = np.linspace(2000, 9000, num=15)
lam2 = h*c / np.linspace(0.1, 10.5, num=105) / 1.6e-19 * 1e9

plt.plot(lam2, data[:,4])
plt.plot(lam2, data[:,6])
plt.plot(lam2, data[:,8])
plt.xlim([min(lam2), 400])
plt.ylim([-23,-16])
plt.show()