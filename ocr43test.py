from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import h5py

f = h5py.File("/SCIA/SDMF31/pieter/vardark_long.h5","r")
ds_vardark = f["varDark"]
ds_orbits = f["dim_orbit"]
idx = (ds_orbits[:] >= 43300) & (ds_orbits[:] <= 43400)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
vd = ds_vardark[idx,0,:]
print(vd.shape)

plt.cla()
imgplot = plt.imshow(vd)
imgplot.set_interpolation('nearest')
plt.show()


f = h5py.File("/SCIA/SDMF31/pieter/interpolated_monthlies_long.h5","r")
ds_aos = f["aos"]
ds_orbits = f["orbits"]
idx = (ds_orbits[:] >= 43300) & (ds_orbits[:] <= 43400)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
aos = ds_aos[idx,:]
print(aos.shape)

plt.cla()
imgplot = plt.imshow(aos)
imgplot.set_interpolation('nearest')
plt.show()

