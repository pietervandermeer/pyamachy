from __future__ import print_function, division

import numpy as np
import h5py
from sciamachy_module import orbitfilter
from vardark_module import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
from scipy.interpolate import interp1d

n_pix = 1024
pixnr = 597
first_orbit = 5000
last_orbit = 6000
#last_orbit = 55000
ofilt = orbitfilter()
ad = AllDarks([1.0, 0.5])

# compute amount of monthlies
n_monthlies = 0
old_monthly = ofilt.get_next_monthly(first_orbit)
while True:
    monthly = ofilt.get_next_monthly(old_monthly)
    if (monthly == old_monthly) or (monthly > last_orbit):
        break
    old_monthly = monthly
    n_monthlies += 1

f = h5py.File("monthly_fits.h5", "w")
orblist = f.create_dataset("orbits", (n_monthlies,), dtype='i')
aos_dset = f.create_dataset("aos", (n_monthlies,n_pix), dtype='f')
lcs_dset = f.create_dataset("lcs", (n_monthlies,n_pix), dtype='f')
amps_dset = f.create_dataset("amps", (n_monthlies,n_pix), dtype='f')
phase_dset = f.create_dataset("channel_phase", (n_monthlies,2), dtype = 'f')
amp2_dset = f.create_dataset("amp2", (n_monthlies,), dtype='f')

old_monthly = ofilt.get_next_monthly(first_orbit)
i_monthly = 0
while True:
    monthly = ofilt.get_next_monthly(old_monthly)
    if (monthly == old_monthly) or (monthly > last_orbit):
        break
    # compute monthly fit, to obtain analog offset
    print("compute montly fit:", monthly, "..")
    channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends = fit_monthly(ad, monthly)
    aos_dset[i_monthly,:] = aos
    lcs_dset[i_monthly,:] = lcs
    amps_dset[i_monthly,:] = amps
    phase_dset[i_monthly,:] = np.array([channel_phase1, channel_phase2])
    amp2_dset[i_monthly] = channel_amp2
    old_monthly = monthly
    orblist[i_monthly] = monthly
    i_monthly += 1

print("done.")

f.close()

