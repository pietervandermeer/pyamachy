import numpy as np
import matplotlib.pyplot as plt
import h5py 

pixnr=571

f30 = h5py.File("/SCIA/SDMF30/sdmf_simudark.h5", "r")
ao30_dset = f30["ch8/ao"]
dc30_dset = f30["ch8/lc"]
orb30_dset = f30["ch8/orbitList"]
ao30_evo = ao30_dset[pixnr,:]
dc30_evo = dc30_dset[pixnr,:]
orbits30 = orb30_dset[:]
f30.close()

#f = h5py.File("/SCIA/SDMF31/pieter/interpolated_monthlies_long_mar2015.h5", "r") 
#ao_dset = f["aos"]
#ao_evo = ao_dset[:,pixnr]
#orbits_dset = f["orbits"]
#orbits = orbits_dset[:]
#f.close()

f = h5py.File("/SCIA/SDMF31/pieter/vardark_long_mar2015.h5", "r") 
dc_dset = f["varDark"]
dc_evo = dc_dset[:,20,pixnr]
orbitsdc_dset = f["dim_orbit"]
orbitsdc = orbitsdc_dset[:]
f.close()

plt.cla()

#dc
plt.title("dark current (BU/s), pixel "+str(pixnr))
plt.xlabel("absolute orbit")
plt.ylabel("dark current (BU/s)")
plt.plot(orbitsdc, dc_evo, 'b.')
plt.plot(orbits30, dc30_evo, 'g.')

plt.show()

