from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py

start = 43300-10
end = 43400+10
pixnr=535

#f = h5py.File("/SCIA/SDMF31/pieter/vardark_long.h5","r")
f = h5py.File("vardark_long.h5","r")
ds_vardark = f["varDark"]
ds_orbits = f["dim_orbit"]
idx = (ds_orbits[:] >= start) & (ds_orbits[:] <= end)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
vd = ds_vardark[idx,0,:]
print(vd.shape)

print(ds_vardark[idx,0,pixnr])

plt.cla()
imgplot = plt.imshow(vd, cmap = cm.Greys_r)
imgplot.set_interpolation('nearest')
plt.show()


#f = h5py.File("/SCIA/SDMF31/pieter/interpolated_monthlies_long.h5","r")
f = h5py.File("interpolated_monthlies_long.h5","r")
ds_aos = f["aos"]
ds_orbits = f["orbits"]
idx = (ds_orbits[:] >= start) & (ds_orbits[:] <= end)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
aos = ds_aos[idx,:]
print(aos.shape)

plt.cla()
imgplot = plt.imshow(aos, cmap = cm.Greys_r)
imgplot.set_interpolation('nearest')
plt.show()

