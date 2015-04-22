import h5py
import matplotlib.pyplot as plt
import numpy as np

import distinct_colours
import envisat

cols = distinct_colours.get_distinct(8)

f = h5py.File("/SCIA/SDMF30/sdmf_pixelmask.h5")
ds_combined = f["smoothMask/combined"]
ds_noise = f["smoothMask/noise"]
ds_invalid = f["smoothMask/invalid"]
ds_residual = f["smoothMask/residual"]
ds_dark = f["smoothMask/darkCurrent"]
ds_chi = f["smoothMask/chiSquare"]
ds_ao = f["smoothMask/analogOffset"]
ds_darkerror = f["smoothMask/darkCurrentError"]
ds_aoerror = f["smoothMask/analogOffsetError"]
ds_orbits = f["smoothMask/orbitList"]

idx = np.argsort(ds_orbits[:])
orbits = ds_orbits[:][idx]
combined = ds_combined[:,:][7*1024:,idx]
invalid = ds_invalid[:,:][7*1024:,idx]
noise = ds_noise[:,:][7*1024:,idx]
residual = ds_residual[:,:][7*1024:,idx]
dark = ds_dark[:,:][7*1024:,idx]
chi = ds_chi[:,:][7*1024:,idx]
darkerror = ds_darkerror[:,:][7*1024:,idx]
aoerror = ds_aoerror[:,:][7*1024:,idx]

combined_sums = combined.sum(axis=0)
noise_sums = noise.sum(axis=0)
invalid_sums = invalid.sum(axis=0)
residual_sums = residual.sum(axis=0)
dark_sums = dark.sum(axis=0)
chi_sums = chi.sum(axis=0)
darkerror_sums = darkerror.sum(axis=0)
aoerror_sums = aoerror.sum(axis=0)

plt.cla()

ax = plt.gca()
ax.set_title("Evolution of SDMF 3.0 pixel mask, channel 8\n\n")
ax.set_xlabel("Orbit number")
ax.set_ylabel("Total bad pixels")

plt.plot(orbits, combined_sums, color=cols[0], label='Total bad pixels')
plt.plot(orbits, invalid_sums, color=cols[1], label='invalid')
plt.plot(orbits, noise_sums, color=cols[2], label='noise')
plt.plot(orbits, residual_sums, color=cols[3], label='residual')
plt.plot(orbits, dark_sums, color=cols[4], label='dark')
plt.plot(orbits, chi_sums, color=cols[5], label='chi^2')
plt.plot(orbits, darkerror_sums, color=cols[6], label='dark error')
plt.plot(orbits, aoerror_sums, color=cols[7], label='analog offset error')

plt.legend(loc='upper left', scatterpoints=10)

ax2 = ax.twiny()
ax2.set_xlabel("Date")
dates = envisat.convert_orbit_to_jd(orbits)
ax2.plot_date(dates,combined_sums,visible=False)

plt.show()


