from __future__ import print_function, division

import numpy as np
import h5py
#from sciamachy_module import orbitfilter
from vardark_module import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
from scipy.interpolate import interp1d

n_pix = 1024
pixnr = 597 # test pixel
# make a graph of pixel "pixnr" instead of computing all pixels
makePiccy = False
# use fitting procedure for every orbit y/n. if no, just find a good trending point by averaging local darks (much faster).
useTrendFit = False
ad = AllDarks([1.0, 0.5])

# open interpolated monthlies
fin = h5py.File("interpolated_monthlies.h5", "r")
orblist = fin["orbits"]
last_orbit = np.max(orblist)
first_orbit = np.min(orblist)
orbit_range = [first_orbit, last_orbit]
inter_aos = fin["aos"]
inter_amps = fin["amps"]
inter_phases = fin["phases"]
inter_amp2 = fin["amp2"]

# get the raw darks
print("get darks..")
n_darks, dummy, pets, coadds, readouts, noise, ephases = ad.get_range(orbit_range)
print("done.")

# sort darks
print("sort darks..")
idx = np.argsort(ephases)
ephases = ephases[idx]
pets = pets[idx]
coadds = coadds[idx]
readouts = readouts[idx,:]
noise = noise[idx,:]
x = ephases
xmin = np.min(x)
xmax = np.max(x)
print("done.")

# subtract interpolated analog offset to get thermal signal, and normalize by time
print("get thermal background signal..")
thermal_background = readouts
i_orbit = 0
m = 0
eorbits = ephases.astype(np.int32)
for orbit in orblist[:]:
    aos = inter_aos[i_orbit,:] 
    n = np.sum(eorbits == orbit)
    if n > 0:
        thermal_background[m:m+n,:] -= (np.matrix(aos).T * np.ones(n)).T
    i_orbit += 1
    m += n
thermal_background /= np.matrix(pets).T * np.ones(n_pix)
y = thermal_background[:,pixnr]
print("done.")

n_tpts = int(xmax) - int(xmin)
trending_phis = np.empty(n_tpts)
trending_ys = np.empty([n_tpts, n_pix])

# determine trending point for each orbit
print("compute trending points..")
i_trend = 0
i_orbit = 0
avg_phi = 0.
xt = np.array([trending_phase]) # trending_phase
lcs = np.ones(n_pix) * 5000 # initial guess for thermal background signal. 5000 BU/s is a good average
for orbit in range(int(xmin), int(xmax)):
    #print(orbit)
    if useTrendFit:
        aos = inter_aos[i_orbit,:]
        amps = inter_amps[i_orbit,:]
        channel_phase1 = inter_phases[i_orbit,0]
        channel_phase2 = inter_phases[i_orbit,1]
        channel_amp2 = inter_amp2[i_orbit,:]
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
    i_orbit += 1

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
xnewi = xnew.astype(np.int32)
print("generating interpolators..")
if makePiccy:
    f = interp1d(trending_phis, trending_ys[:,pixnr])
    f2 = interp1d(trending_phis, trending_ys[:,pixnr], kind='cubic')
else:
    print(trending_phis.shape, trending_ys.shape)
    f = interp1d(trending_phis, trending_ys, axis=0)
print("done.")
print("interpolate trend..")
xt = np.array([phi_t]) # trending_phase
wave = np.empty([xnew.size, n_pix])
out_orblist = np.array([], dtype=np.int32)
for orbit in range(int(ximin), int(ximax)):
    out_orblist = np.append(out_orblist, orbit)
    i_orbit = (np.where(orblist[:] == orbit))[0][0]
    print(i_orbit, orblist[0], orbit)
    aos = inter_aos[i_orbit,:]
    amps = inter_amps[i_orbit,:]
    channel_phase1 = inter_phases[i_orbit,0]
    channel_phase2 = inter_phases[i_orbit,1]
    channel_amp2 = inter_amp2[i_orbit]
    idx = xnewi == orbit
    xnew_ = xnew[idx]
    if makePiccy:
        p = aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phase1, channel_amp2, channel_phase2
        wave_ = scia_dark_fun2n(p, xnew_) - scia_dark_fun2n(p, xt)
        wave[idx, pixnr] = wave_ 
    else:
        for i_pix in range(n_pix):
            p = aos[i_pix], lcs[i_pix], amps[i_pix], 0, channel_phase1, channel_amp2, channel_phase2
            wave_ = scia_dark_fun2n(p, xnew_) - scia_dark_fun2n(p, xt)
            wave[idx, i_pix] = wave_ 
fin.close() # close interpolated monthlies databse, won't be needing it anymore
print("done.")

if makePiccy:
    # single pixel
    wave[:,pixnr] += f(xnew)
    import matplotlib.pyplot as plt
    plt.ticklabel_format(useOffset=False)
    plt.plot(x,y,'v', trending_phis, trending_ys[:,pixnr], 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--', xnew, wave[:,pixnr],'-')
    #plt.plot(trending_phis, trending_ys, 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
    plt.legend(['orig data', 'avg data', 'linear', 'cubic', 'reconstruct'], loc='best')
    plt.show()
else:
    # store all pixels in output hdf5
    print("storing..")
    wave += f(xnew)
    fout = h5py.File("vardark_long.h5", "w")
    print(out_orblist.size)
    orb_dset = fout.create_dataset("orbits", out_orblist.shape, dtype='i')
    orb_dset[:] = out_orblist
    wave_dset = fout.create_dataset("wave", wave.shape, dtype='f')
    wave_dset[:,:] = wave
    fout.close()    
    print("done.")

