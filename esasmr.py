from __future__ import print_function, division

import h5py
import numpy as np
import matplotlib.pyplot as plt

#orbit = 52857
#orbit = 14000
orbit = 52000

#
# ESA
#

f = h5py.File("/SCIA/SDMF31/sdmf_extract_l1b.h5","r")
ds_sun = f["W/sun_diffuser"]
ds_orbits = f["W/orbitList"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("esa idx", idx)
esa_smr = (ds_sun[idx])[7*1024:]

#
# SDMF SMR 3.2
#

#f = h5py.File("/SCIA/SDMF31/pieter/sdmf_smr_123_.h5","r")
#f = h5py.File("/SCIA/SDMF31/pieter/sdmf_smr_raw_.h5","r")
f = h5py.File("/SCIA/SDMF31/pieter/sdmf_smr_12_.h5","r")
#f = h5py.File("/SCIA/SDMF31/pieter/sdmf_smr_full.h5","r")
#f = h5py.File("smr_14009_1237.h5","r")
#f = h5py.File("smr_14009_12345678.h5","r")
ds_orbits = f["orbitList"]
ds_sun = f["smr"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
smr32_orbit = ds_orbits[idx]
print("sdmf32 idx", idx, "smr orbit", smr32_orbit)
sdmf_smr = (ds_sun[idx])[7*1024:] # * 5e9

#
# SDMF SMR 3.1
#

f = h5py.File("/SCIA/SDMF31/sdmf_smr.h5","r")
ds_orbits = f["orbitList"]
ds_sun = f["smr"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
smr31_orbit = ds_orbits[idx]
print("sdmf31 idx", idx, "smr orbit", smr31_orbit)
sdmf31_smr = (ds_sun[idx])[7*1024:] # * 5e9

#
# SDMF SMR 3.0
#

f = h5py.File("/SCIA/SDMF30/sdmf_smr.h5","r")
ds_orbits = f["orbitList"]
ds_sun = f["SMR"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
smr30_orbit = ds_orbits[idx]
print("sdmf30 idx", idx, "smr orbit", smr30_orbit)
sdmf30_smr = ds_sun[7*1024:,idx] # * 5e9

f = h5py.File("/SCIA/SDMF31/pieter/sdmf_smooth_pyxelmask32_.h5")
ds_orbits = f["orbits"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("mask idx", idx)
bad_mask = f["combinedFlag"][idx,:]
bad_mask = np.logical_or(bad_mask, (sdmf_smr < 100)) # also kick out all the very low smr pixels
good_idx = np.where(np.logical_not(bad_mask))
print(good_idx[0].shape, sdmf_smr.shape)
good_idx = good_idx[0]
print(good_idx, good_idx.size)  
all_idx = np.arange(1024)

#chosen_idx = all_idx
chosen_idx = good_idx

#for val in sdmf_smr[all_idx]:
#	print(val)

plt.cla()
plt.ticklabel_format(useOffset=False)
plt.title("orbit "+str(orbit))
#plt.ylim([0,10000000000])
#plt.plot(chosen_idx, esa_smr[chosen_idx], 'b.', label="ESA")
plt.plot(all_idx, sdmf_smr[all_idx], 'g.', label="SDMF3.2 all")
plt.plot(good_idx, sdmf_smr[good_idx], 'r.', label="SDMF3.2 good")
plt.plot(good_idx, sdmf30_smr[good_idx], 'k.', label="SDMF3.0 good")
plt.plot(good_idx, sdmf31_smr[good_idx], 'c.', label="SDMF3.1 good")
plt.legend(loc="best")
plt.show()

