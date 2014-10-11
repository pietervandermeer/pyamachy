from __future__ import print_function, division

import numpy as np
from sciamachy_module import orbitfilter
from vardark_module import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
from scipy.interpolate import interp1d

# use fitting procedure for every orbit y/n. if no, just find a good trending point by averaging local darks (much faster).
useTrendFit = False
n_pix = 1024
ofilt = orbitfilter()
#orbit_range = [26793, 27232]
orbit_range = [27228, 27745]
#orbit_range = [27016, 27030]
monthly = ofilt.get_next_monthly(orbit_range[0])
pixnr = 597

# import darks
print("importing darks..")
ad = AllDarks([1.0, 0.5])
print("done.")

# get the raw darks
print("get darks..")
n_darks, dummy, pets, coadds, readouts, noise, ephases = ad.get_range(orbit_range)
print("done.")

# compute monthly fit, to obtain analog offset
print("compute fit..")
channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends = fit_monthly(ad, monthly)
#amps*=1.1 # hehehe
print("done.")

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

# subtract analog offset to get thermal signal, and normalize by time
thermal_background = readouts
thermal_background -= (np.matrix(aos).T * np.ones(n_darks)).T
thermal_background /= np.matrix(pets).T * np.ones(n_pix)
y = thermal_background[:,pixnr]

n_tpts = int(xmax) - int(xmin)
trending_phis = np.empty(n_tpts)
trending_ys = np.empty([n_tpts, n_pix])

# determine trending point for each orbit
print("compute trending points..")
i_trend = 0
avg_phi = 0.
xt = np.array([trending_phase]) # trending_phase
for orbit in range(int(xmin), int(xmax)):
    #print(orbit)
    if useTrendFit:
        x__, lcs_, res_trends, dum1, dum2 = fit_eclipse_orbit(ad, orbit, aos, lcs, amps, channel_amp2, channel_phase1, channel_phase2)
        for i_pix in range(n_pix):
            p = aos[i_pix], lcs_[i_pix], amps[i_pix], res_trends[i_pix], channel_phase1, channel_amp2, channel_phase2
            trending_ys[i_trend, i_pix] = scia_dark_fun2m(p, xt)
        avg_phi += trending_phase
        trending_phis[i_trend] = trending_phase+orbit
        print(aos[pixnr], lcs_[pixnr], amps[pixnr], res_trends[pixnr], channel_phase1, channel_amp2, channel_phase2)
        print(trending_ys[i_trend, pixnr])
        i_trend += 1
    else:
        idx = (ephases >= orbit) & (ephases < (orbit+1))
        local_phi = ephases[idx]
        local_y = thermal_background[idx,:]
        #print(local_phi)
        abs_dist1 = np.abs(local_phi - (trending_phase+orbit))
        idx1 = (np.argsort(abs_dist1))[0:6]
        phi1 = np.mean(local_phi[idx1])
        #print(local_y.shape)
        avg1 = np.mean(local_y[idx1,:], axis=0)
        trending_phis[i_trend] = phi1
        trending_ys[i_trend,:] = avg1
        if not np.isnan(phi1):
            avg_phi += (phi1%1)
            i_trend += 1
        #print(phi1%1)

avg_phi /= i_trend
print(avg_phi, i_trend)
trending_phis = trending_phis[0:i_trend]
trending_ys = trending_ys[0:i_trend]
print("done.")

print("remove invalid entries..")
idx_goodphi = np.isfinite(trending_phis)
trending_phis = trending_phis[idx_goodphi]
trending_ys = trending_ys[idx_goodphi,:]
# filter a single pixel. probably a bad idea
#idx_goody = np.isfinite(trending_ys[:,pixnr])
#trending_phis = trending_phis[idx_goody]
#trending_ys = trending_ys[idx_goody,:]
print("done.")

ximin = int(xmin)+.9
ximax = int(xmax)-.9
xnew = np.linspace(ximin, ximax, (ximax-ximin)*50.)
phi_t = avg_phi # trending_phase
print("phi_t=",phi_t)
#print(trending_ys[:,pixnr], trending_phis)

# visualize a pixel of choice
#for ytje in trending_ys[:,pixnr]:
#    print(ytje)
#for phitje in trending_phis:
#    print(phitje)
f = interp1d(trending_phis, trending_ys[:,pixnr])
f2 = interp1d(trending_phis, trending_ys[:,pixnr], kind='cubic')
p = aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phase1, channel_amp2, channel_phase2
xt = np.array([phi_t]) # trending_phase
#wave = scia_dark_fun2n(p, xnew) + f2(xnew) - scia_dark_fun2n(p, xt) 
wave = scia_dark_fun2n(p, xnew) + f(xnew) - scia_dark_fun2n(p, xt) 
print(xnew, wave, scia_dark_fun2n(p, xnew), scia_dark_fun2n(p, xt), f(xnew))

import matplotlib.pyplot as plt
plt.ticklabel_format(useOffset=False)
plt.plot(x,y,'v', trending_phis, trending_ys[:,pixnr], 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--', xnew, wave,'-')
#plt.plot(trending_phis, trending_ys, 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
plt.legend(['orig data', 'avg data', 'linear', 'cubic', 'reconstruct'], loc='best')
plt.show()

print("interpolate trend..")
xt = np.array(phi_t) # trending_phase
for i_pix in range(n_pix):
    print(i_pix)
    #f2 = interp1d(trending_phis, trending_ys[:,i_pix], kind='cubic')
    f = interp1d(trending_phis, trending_ys[:,i_pix])
    p = aos[i_pix], lcs[i_pix], amps[i_pix], 0, channel_phase1, channel_amp2, channel_phase2
    #wave = scia_dark_fun2n(p, xnew) + f2(xnew) - scia_dark_fun2n(p, xt)
    wave = scia_dark_fun2n(p, xnew) + f(xnew) - scia_dark_fun2n(p, xt)
print("done.")


