import h5py
import numpy as np
import matplotlib.pyplot as plt

orbit = 14009

f = h5py.File("/SCIA/SDMF31/sdmf_extract_l1b.h5","r")
ds_sun = f["W/sun_diffuser"]
ds_orbits = f["W/orbitList"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("esa idx", idx)
esa_smr = (ds_sun[idx])[7*1024:]

#f = h5py.File("sdmf_smr.h5","r")
f = h5py.File("smr_test.h5","r")
ds_orbits = f["orbitList"]
ds_sun = f["smr"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("sdmf idx", idx)
sdmf_smr = (ds_sun[idx])[7*1024:] # * 5e9

f = h5py.File("sdmf_smooth_pyxelmask32.h5")
ds_orbits = f["orbits"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("mask idx", idx)
good_idx = np.where(np.logical_not(f["combinedFlag"][idx,:]))
good_idx = good_idx[0]
print good_idx, good_idx.size  
all_idx = np.arange(1024)

#chosen_idx = all_idx
chosen_idx = good_idx

plt.cla()
plt.ticklabel_format(useOffset=False)
plt.plot(chosen_idx, sdmf_smr[chosen_idx], 'r.', label="SDMF3.2")
plt.plot(chosen_idx, esa_smr[chosen_idx], 'b.', label="ESA")
plt.legend(loc="best")
plt.show()
