from __future__ import print_function, division

import numpy as np
import h5py
import matplotlib.pyplot as plt

fname = "/SCIA/SDMF31/sdmf_extract_calib.h5"
f = h5py.File(fname, "r")

orbits = f["/State_08/orbitList"]
mtbl = f["/State_08/metaTable"]
plot_x = orbits[:]
tdets = mtbl[:]['detectorTemp']
tobm = mtbl[:]['obmTemp']
print(tdets.shape)
plot_y = tdets[:,7]
plot_y2 = tobm[:]

plt.ticklabel_format(useOffset=False)
print(plot_x.shape, plot_y.shape)
plt.plot(plot_x,plot_y,'o',plot_x,plot_y2,'v')
plt.show()
