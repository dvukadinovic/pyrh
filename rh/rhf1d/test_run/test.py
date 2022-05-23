import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import time
import copy
import sys

import globin
import pyrh

atmos = globin.Atmosphere("atmos_combined_ss_v2.fits")

start = time.time()

aux = pyrh.RH()
# aux.dummy()
# sys.exit()

idx, idy = 0, 0
ne, nH = aux.hse(0, atmos.data[idx, idy, 0], atmos.data[idx, idy, 1], atmos.data[idx, idy, 2],
                   atmos.data[idx, idy, 3], atmos.data[idx, idy, 4],
                   atmos.data[idx, idy, 5] / 1e4, atmos.data[idx, idy, 6], atmos.data[idx, idy, 7],
                   atmos.data[idx, idy, 8:], 0, fudge_lam, fudge)
atmos.data[idx,idy,2] = ne/1e6
atmos.data[idx,idy,8:] = nH/1e6
# plt.plot(atmos.data[idx,idy,9])
# plt.plot(nH[1]/1e6)
# plt.yscale("log")
# plt.show()
# print(time.time() - start)

# sys.exit()

loggf_ids = np.array([], dtype=np.int32)
loggf_values = np.array([], dtype=np.float64)

lam_ids = np.array([], dtype=np.int32)
lam_values = np.array([], dtype=np.float64)

# def spec(args):
idx, idy = 0, 0
spec = aux.compute1d(0, atmos.data[idx, idy, 0], atmos.data[idx, idy, 1], atmos.data[idx, idy, 2],
                   atmos.data[idx, idy, 3], atmos.data[idx, idy, 4],
                   atmos.data[idx, idy, 5] / 1e4, atmos.data[idx, idy, 6], atmos.data[idx, idy, 7],
                   atmos.data[idx, idy, 8:], 1, fudge_lam, fudge,
                   loggf_ids, loggf_values,
                   lam_ids, lam_values/1e4)

# pool = mp.Pool(4)
#
# idx = np.arange(atmos.nx)
# idy = np.arange(atmos.ny)
# IDx, IDy = np.meshgrid(idx, idy)
# args = zip(IDx.flatten(), IDy.flatten())
# pool.map(spec, iterable=args)

print(time.time() - start)
obs = globin.Observation("obs.fits")

plt.plot(spec.I)
plt.plot(obs.spec[0,0,:,0])

# plt.plot(obs.spec[0,0,:,0] - spec.I)

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
