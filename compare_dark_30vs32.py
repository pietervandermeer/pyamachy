from __future__ import print_function, division

import sys
import matplotlib.pyplot as plt
import numpy as np 
import h5py

#orbit = 9000
#orbit = 51990
orbit = 52300

#
# read SDMF 3.2 data
#

#vd_fname = "/SCIA/SDMF31/pieter/vardark_long__.h5"
#ao_fname = "/SCIA/SDMF31/pieter/interpolated_monthlies_long__.h5"
vd_fname = "/SCIA/SDMF31/pieter/vardark_long_mar2015.h5"
ao_fname = "/SCIA/SDMF31/pieter/interpolated_monthlies_long_mar2015.h5"
vd_fid = h5py.File(vd_fname,"r")
ao_fid = h5py.File(ao_fname,"r")
vd_dset = vd_fid["varDark"]
vd_orb_idx = vd_fid["dim_orbit"][:] == orbit
if (np.any(vd_orb_idx)):
	vd_orb_idx = np.where(vd_orb_idx)[0]
	print(vd_orb_idx)
	#vd = vd_dset[vd_orb_idx, 10, :].flatten()
	sl = vd_dset[vd_orb_idx, 10, :]
	vd = sl.flatten()
else:
	print("orbit "+str(orbit)+" not found in "+vd_fname)
	sys.exit(1)

ao_dset = ao_fid["aos"]
ao_orb_idx = ao_fid["orbits"][:] == orbit
if (np.any(ao_orb_idx)):
	ao_orb_idx = np.where(ao_orb_idx)[0]
	print(ao_orb_idx)
	ao = (ao_dset[ao_orb_idx, :]).flatten()
else:
	print("orbit "+str(orbit)+" not found in "+ao_fname)
	sys.exit(1)

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
plt.plot(pix_axis, ao+vd*.5, 'b.', label="3.2")
plt.plot(pix_axis, dark, 'r.', label="3.0")
plt.legend()
plt.show()
