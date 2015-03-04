from __future__ import print_function, division

import matplotlib.pyplot as plt
import numpy as np 
import h5py

orbit = 9000

#
# read SDMF 3.2 data
#

vd_fname = "/SCIA/SDMF31/pieter/vardark_long__.h5"
ao_fname = "/SCIA/SDMF31/pieter/interpolated_monthlies_long__.h5"
vd_fid = h5py.File(vd_fname,"r")
ao_fid = h5py.File(ao_fname,"r")
vd_dset = vd_fid["varDark"]
vd_orb_idx = vd_fid["dim_orbit"][:] == orbit
vd = vd_dset[vd_orb_idx, 10, :].flatten()
ao_dset = ao_fid["aos"]
ao_orb_idx = ao_fid["orbits"][:] == orbit
ao = ao_dset[ao_orb_idx, :].flatten()
print(ao.shape,vd.shape)

#
# read SDMF 3.0 data
#

dark_fname = "/SCIA/SDMF30/sdmf_simudark.h5"
dark_fid = h5py.File(dark_fname,"r")
ao_dset = dark_fid["ch8/ao"]
lc_dset = dark_fid["ch8/lc"]
orb30_idx = dark_fid["ch8/orbitList"][:] == orbit
dark = (ao_dset[:,orb30_idx] + 0.5*lc_dset[:,orb30_idx]).flatten() # no orbvar yet

#
# plot 
#

pix_axis = np.arange(1024)
plt.cla()
plt.plot(pix_axis, ao+vd*.5, 'b.')
plt.plot(pix_axis, dark, 'g.')
plt.show()
