#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt 
import numpy as np 
import h5py 

start_orbit = 47310
end_orbit = 47400

#f = h5py.File("sdmf_smooth_pyxelmask32.h5", "r")
f = h5py.File("/SCIA/SDMF31/pieter/sdmf_pyxelmask32.h5", "r")
ds_orbits = f["orbits"]
ds_sat = f["saturation"]
idx = np.argmin(np.abs(ds_orbits[:] - start_orbit))
orb_data = ds_sat[idx]

plt.cla()
plt.plot(orb_data, 'b.')
plt.show()
