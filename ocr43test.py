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
ds_uncertainty = f["uncertainties"]
ds_orbits = f["dim_orbit"]
idx = (ds_orbits[:] >= start) & (ds_orbits[:] <= end)
if np.sum(idx) > 0:
    orbits_vd = ds_orbits[idx]
    idx = np.where(idx)[0]
else:
    raise Exception("No orbits found")
vd = ds_vardark[idx,phasebin,:]
un = ds_uncertainty[idx,:]
print(vd.shape)

print(ds_vardark[idx,phasebin,pixnr])

# plt.cla()
# imgplot = plt.imshow(vd, cmap = cm.Greys_r)
# imgplot.set_interpolation('nearest')
# plt.show()


f = h5py.File("/SCIA/SDMF31/pieter/interpolated_monthlies_long_.h5","r")
#f = h5py.File("interpolated_monthlies_long.h5","r")
ds_aos = f["aos"]
ds_orbits = f["orbits"]
idx = (ds_orbits[:] >= start) & (ds_orbits[:] <= end)
if np.sum(idx) > 0:
    idx = np.where(idx)[0]
else:
    raise Exception("No orbits found")
aos = ds_aos[idx,:]
print(aos.shape)

if False:
    #2d plots
    plt.cla()
    fig = plt.figure()
    ax = fig.add_subplot(3,1,1)
    imgplot1 = ax.imshow(np.log(vd), cmap=cm.Greys_r, aspect="auto")
    imgplot1.set_interpolation('nearest')
    ax2 = fig.add_subplot(3,1,2, sharex=ax, sharey=ax)
    imgplot2 = ax2.imshow(np.log(un), cmap=cm.Greys_r, aspect="auto")
    imgplot2.set_interpolation('nearest')
    ax3 = fig.add_subplot(3,1,3, sharex=ax, sharey=ax)
    imgplot3 = ax3.imshow(np.log(aos), cmap=cm.Greys_r, aspect="auto")
    imgplot3.set_interpolation('nearest')
    plt.show()

#1d
pixnr = 515
norm_dark = vd[:,pixnr]
norm_dark -= np.mean(norm_dark)
norm_dark /= np.max(np.abs(norm_dark))
norm_un = un[:,pixnr]
idx = np.logical_not(np.isfinite(norm_un))
print(np.where(idx))
norm_un = np.nan_to_num(norm_un)
norm_un -= np.mean(norm_un)
norm_un /= np.max(np.abs(norm_un))
quo = vd[:,pixnr]/un[:,pixnr]

plt.cla()
plt.ticklabel_format(useOffset=False)
plt.plot(orbits_vd, norm_dark, 'b-')
plt.plot(orbits_vd, norm_un, 'r-')
plt.plot(orbits_vd, quo, 'g-')
plt.show()

