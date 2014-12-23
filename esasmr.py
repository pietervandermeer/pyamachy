import h5py
import numpy as np
import matplotlib.pyplot as plt

orbit = 43000

f = h5py.File("/SCIA/SDMF31/sdmf_extract_l1b.h5","r")
ds_sun = f["W/sun_diffuser"]
ds_orbits = f["W/orbitList"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("esa idx", idx)
esa_smr = (ds_sun[idx])[7*1024:]

f = h5py.File("sdmf_smr_123.h5","r")
ds_orbits = f["orbitList"]
ds_sun = f["smr"]
idx = np.argmin(np.abs(ds_orbits[:] - orbit))
print("sdmf idx", idx)
sdmf_smr = (ds_sun[idx])[7*1024:] * 5e9

plt.cla()
plt.ticklabel_format(useOffset=False)
plt.plot(esa_smr, 'bo', label="ESA")
plt.plot(sdmf_smr, 'ro', label="SDMF3.2")
plt.legend(loc="best")
plt.show()
