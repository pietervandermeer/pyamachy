from __future__ import print_function, division

import h5py
import numpy as np

if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    fid = h5py.File("/SCIA/SDMF30/sdmf_pixelmask.h5")
    dset_orbit = fid["orbitalMask/orbitList"]
    dset_wls = fid["orbitalMask/wlsResponse"]
    dset_sun = fid["orbitalMask/sunResponse"]
    dset_ppg = fid["orbitalMask/pixelGain"]

    start_orbit = 45000
    end_orbit = 50000

    orbit_idx = (dset_orbit[:] >= start_orbit) & (dset_orbit[:] <= end_orbit)

    print(dset_wls.shape)
    print(dset_sun.shape)
    print(dset_ppg.shape)

    orbits = dset_orbit[orbit_idx]
    i = 0
    for orbit in orbits:
#    for orbit in range(0,1000):
        #msk = dset_sun[7*1024:,i]
        #msk = dset_wls[7*1024:,i]
        msk = dset_ppg[7*1024:,i]
        print(orbit, np.where(msk))
        #print(orbit, dset[7*1024:,orbit])
        i += 1
