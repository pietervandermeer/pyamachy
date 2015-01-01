import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File("sdmf_smooth_pyxelmask32.h5","r")
plt.cla()
ds_orbits = f["orbits"]
ds_combf = f["combinedFlag"]
ds_err = f["darkError"]
idx = np.argmin(np.abs(ds_orbits[:] - 43300))
mask1 = ds_combf[idx,:]
err1 = ds_err[idx,:]
idx = np.argmin(np.abs(ds_orbits[:] - 43400))
mask2 = ds_combf[idx,:]
err2 = ds_err[idx,:]
plt.cla()
plt.plot(err1,"b.")
plt.plot(err2,"r.")
plt.show()

