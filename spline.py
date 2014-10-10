from __future__ import print_function

import numpy as np
from sciamachy_module import orbitfilter
from vardark_module import AllDarks, trending_phase, fit_monthly
from scipy.interpolate import interp1d

ofilt = orbitfilter()
orbit_range = [27016,27030]
monthly = ofilt.get_next_monthly(orbit_range[0])
pixnr = 597

# import darks
print("importing darks..")
ad = AllDarks([1.0, 0.5])
ad.lump(orbit_range)
ad.finalize()
print("done.")

# get the raw darks
n_darks, dummy, pets, coadds, readouts, noise, ephases = ad.get_range(orbit_range)

# compute monthly fit, to obtain analog offset
channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends = fit_monthly(ad, monthly)

# sort darks
idx = np.argsort(ephases)
ephases = ephases[idx]
pets = pets[idx]
coadds = coadds[idx]
readouts = readouts[idx,:]
noise = noise[idx,:]
x = ephases
xmin = np.min(x)
xmax = np.max(x)

# take one pixel and normalize its readouts
y = readouts[:, pixnr].flatten() - aos[pixnr]
upets = np.unique(pets)
minpet = np.min(upets)
maxpet = np.max(upets)
idx05 = (pets == minpet)
idx10 = (pets == maxpet)
y[idx05] /= minpet
y[idx10] /= maxpet

n_tpts = int(xmax) - int(xmin)
trending_phis = np.empty(n_tpts)
trending_ys = np.empty(n_tpts)

# determine trending point for each orbit
i_trend = 0
for orbit in range(int(xmin), int(xmax)):
    print(orbit)
    idx = (ephases >= orbit) & (ephases < (orbit+1))
    local_phi = ephases[idx]
    local_y = y[idx]
    print(local_phi)
    abs_dist1 = np.abs(local_phi - (trending_phase+orbit))
    idx1 = (np.argsort(abs_dist1))[0:3]
    phi1 = np.mean(local_phi[idx1])
    avg1 = np.mean(local_y[idx1])
    print(phi1, avg1)
    trending_phis[i_trend] = phi1
    trending_ys[i_trend] = avg1
    i_trend += 1

f = interp1d(trending_phis, trending_ys)
f2 = interp1d(trending_phis, trending_ys, kind='cubic')

ximin = int(xmin)+1
ximax = int(xmax)-1
xnew = np.linspace(ximin, ximax, (ximax-ximin)*100.)
#for x_ in x:
#    print(x_)
#for x_ in xnew:
#    print(x_)
import matplotlib.pyplot as plt
plt.plot(x,y,'-', trending_phis, trending_ys, 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
#plt.plot(trending_phis, trending_ys, 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
plt.legend(['orig data', 'avg data', 'linear', 'cubic'], loc='best')
plt.show()
