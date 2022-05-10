import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import time
import sys
import pickle

import globin
import pyrh

atmos = globin.Atmosphere("atmos_combined_ss_v2.fits")

start = time.time()


aux = pyrh.RH()

fudge_num = 5
fudge_lam = np.linspace(401.5, 401.7, num=fudge_num)
fudge = np.ones((3, fudge_num))

# def spec(args):
idx, idy = 0, 0
spec = aux.compute1d(0, atmos.data[idx, idy, 0], atmos.data[idx, idy, 1], atmos.data[idx, idy, 2],
                   atmos.data[idx, idy, 3], atmos.data[idx, idy, 4],
                   atmos.data[idx, idy, 5] / 1e4, atmos.data[idx, idy, 6], atmos.data[idx, idy, 7],
                   atmos.data[idx, idy, 8:], 1, fudge_lam, fudge)

# pool = mp.Pool(4)
#
# idx = np.arange(atmos.nx)
# idy = np.arange(atmos.ny)
# IDx, IDy = np.meshgrid(idx, idy)
# args = zip(IDx.flatten(), IDy.flatten())
# pool.map(spec, iterable=args)

print(time.time() - start)

plt.plot(spec.I[:-1])

spec = globin.Observation("obs.fits")
plt.plot(spec.spec[0,0,:,0])

plt.show()
sys.exit()

start = time.time()

run_name = "dummy"
globin.read_input(run_name=run_name)

if globin.mode == 0:
    spec = globin.make_synthetic_observations(globin.atm, globin.noise)

print(time.time() - start)

plt.plot(spec.spec[0, 0, :, 0])
# plt.plot(spec.spec[0,0,:,0] - specRH.I[:-1])

plt.show()
