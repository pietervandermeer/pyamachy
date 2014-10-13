from __future__ import print_function, division

import numpy as np
import h5py
from sciamachy_module import orbitfilter
from vardark_module import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
from scipy.interpolate import interp1d
from scipy import array

#-- globals --------------------------------------------------------------------

n_pix = 1024
first_orbit = 5000
last_orbit = 6000

#-- functions ------------------------------------------------------------------

# extrapolate outside of the original data range of the specified interpolator
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] - (x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0]) # mirror!
        elif x > xs[-1]:
            return ys[-1] - (x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2]) # mirror!
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike

#-- main -----------------------------------------------------------------------

fin = h5py.File("monthly_fits.h5", "r")
monthlies = fin["orbits"]
#n_orbits = monthlies[-1] - monthlies[0]
#orblist = np.arange(monthlies[0], monthlies[-1])
n_orbits = 53200 - monthlies[0] # also include last orbits, beyond last monthly calibration orbit
orblist = np.arange(monthlies[0], 53200)
aos_dset = fin["aos"]
lcs_dset = fin["lcs"]
amps_dset = fin["amps"]
phase_dset = fin["channel_phase"]
amp2_dset = fin["amp2"]

#extend beyond

fout = h5py.File("interpolated_monthlies.h5", "w")
iorbits = fout.create_dataset("orbits", (n_orbits,), dtype='i')
iphases = fout.create_dataset("phases", (n_orbits,2), dtype='f')
iaos = fout.create_dataset("aos", (n_orbits,n_pix), dtype='f')
iamps = fout.create_dataset("amps", (n_orbits,n_pix), dtype='f')
iamp2 = fout.create_dataset("amp2", (n_orbits,), dtype='f')

print(orblist.shape, iphases.shape)
fun_phase1 = extrap1d(interp1d(monthlies, phase_dset[:,0]))
fun_phase2 = extrap1d(interp1d(monthlies, phase_dset[:,1]))
fun_amp2 = extrap1d(interp1d(monthlies, amp2_dset[:]))
# for i_pix in range(n_pix):
#     print(i_pix)
#     fun_aos = extrap1d(interp1d(monthlies, aos_dset[:,i_pix]))
#     fun_amps = extrap1d(interp1d(monthlies, amps_dset[:,i_pix]))
#     iaos[:,i_pix] = fun_aos(orblist)
#     iamps[:,i_pix] = fun_amps(orblist)
fun_aos = extrap1d(interp1d(monthlies, aos_dset[:,:], axis=0))
fun_amps = extrap1d(interp1d(monthlies, amps_dset[:,:], axis=0))
iaos[:,:] = fun_aos(orblist)
iamps[:,:] = fun_amps(orblist)

iorbits[:] = orblist
iphases[:,0] = fun_phase1(orblist)
iphases[:,1] = fun_phase2(orblist)
iamp2[:] = fun_amp2(orblist)

# write out and close
fout.close()
# close input
fin.close()

