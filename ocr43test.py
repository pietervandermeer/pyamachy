#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Display vardark product (DC and AO) over a large orbit range, using a 2D image.
"""

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py

#start = 4500
#end = 5000
start = 43000
end = 44000
#start = 43300
#end = 43400

pixnr = 535
phasebin = 10

#f = h5py.File("/SCIA/SDMF31/pieter/vardark_short_.h5","r")
f = h5py.File("/SCIA/SDMF31/pieter/vardark_long_.h5","r")
#f = h5py.File("vardark_long.h5","r")
ds_vardark = f["varDark"]
ds_orbits = f["dim_orbit"]
idx = (ds_orbits[:] >= start) & (ds_orbits[:] <= end)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
else:
    raise Exception("No orbits found")
vd = ds_vardark[idx,phasebin,:]
print(vd.shape)

print(ds_vardark[idx,phasebin,pixnr])

# plt.cla()
# imgplot = plt.imshow(vd, cmap = cm.Greys_r)
# imgplot.set_interpolation('nearest')
# plt.show()


#f = h5py.File("/SCIA/SDMF31/pieter/interpolated_monthlies_long.h5","r")
f = h5py.File("interpolated_monthlies_long.h5","r")
ds_aos = f["aos"]
ds_orbits = f["orbits"]
idx = (ds_orbits[:] >= start) & (ds_orbits[:] <= end)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
else:
    raise Exception("No orbits found")
aos = ds_aos[idx,:]
print(aos.shape)

# plt.cla()
# imgplot = plt.imshow(aos, cmap = cm.Greys_r)
# imgplot.set_interpolation('nearest')
# plt.show()

plt.cla()
fig = plt.figure()
ax = fig.add_subplot(2,1,1)
imgplot1 = ax.imshow(vd, cmap=cm.Greys_r, aspect="auto")
imgplot1.set_interpolation('nearest')
ax2 = fig.add_subplot(2,1,2, sharex=ax, sharey=ax)
imgplot2 = ax2.imshow(aos, cmap=cm.Greys_r, aspect="auto")
imgplot2.set_interpolation('nearest')
plt.show()
