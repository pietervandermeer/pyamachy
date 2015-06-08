import numpy as np
import matplotlib.pyplot as plt
import h5py 
import matplotlib.cm as cm

winstart=0
winend=1024
#winstart=380
#winend=621
winsize = winend - winstart

f30 = h5py.File("/SCIA/SDMF30/sdmf_simudark.h5", "r")
ao30_dset = f30["ch8/ao"]
dc30_dset = f30["ch8/lc"]
orb30_dset = f30["ch8/orbitList"]
#ao30_evo = ao30_dset[pixnr,:]
dc30_evo = dc30_dset[winstart:winend,:]
orbits30 = orb30_dset[:]
f30.close()

#f = h5py.File("/SCIA/SDMF31/pieter/interpolated_monthlies_long_mar2015.h5", "r") 
#ao_dset = f["aos"]
#ao_evo = ao_dset[winstart:winend,pixnr]
#orbits_dset = f["orbits"]
#orbits = orbits_dset[:]
#f.close()

f = h5py.File("/SCIA/SDMF31/pieter/vardark_long_mar2015.h5", "r") 
dc_dset = f["varDark"]
dc_evo = dc_dset[:,20,winstart:winend]
orbitsdc_dset = f["dim_orbit"]
orbitsdc = orbitsdc_dset[:]
f.close()

orbit32_idx = np.in1d(orbitsdc,orbits30)
orbit30_idx = np.in1d(orbits30,orbitsdc)
dat32_ = dc_evo[orbit32_idx,:]
dat30_ = dc30_evo[:,orbit30_idx].transpose()
diff_ = np.log(np.abs(dat32_ - dat30_))
diff_aligned = np.zeros((winsize,53200))
i = 0
for orbit in orbitsdc[orbit32_idx]:
	diff_aligned[:,orbit] = diff_[i,:]
	i += 1

print str(orbitsdc[orbit32_idx][19355])
print str(orbitsdc[orbit32_idx][23500])
print str(orbitsdc[orbit32_idx][24000])
print str(orbitsdc[orbit32_idx][25000])

plt.cla()
plt.title("difference between SDMF3.0 and SDMF3.2 dark signal (BU/s)")
plt.xlabel("absolute orbit")
plt.ylabel("dark current difference (BU/s)")
plt.imshow(diff_aligned, cmap = cm.Greys_r, interpolation='none')
plt.show()

