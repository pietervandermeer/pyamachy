from __future__ import print_function, division

import numpy as np 
import h5py
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.misc
import warnings
warnings.simplefilter("error") # warnings to errors

# smoothing window size
winsize = 50

#-- functions ------------------------------------------------------------------

def smooth(orbits, data):
    smoothed = np.empty(data.shape, dtype=np.float64)
    num_orbits = orbits.size
    for i_orbit in range(num_orbits):
        if i_orbit < winsize/2:
            window = np.arange(0, i_orbit+winsize/2, dtype=np.int)
        elif i_orbit >= num_orbits - winsize/2:
            window = np.arange(i_orbit-winsize/2, num_orbits, dtype=np.int)
        else:
            window = np.arange(i_orbit-winsize/2, i_orbit+winsize/2, dtype=np.int)
        smoothed[i_orbit, :] = np.mean(data[window,:], axis=0)
    return smoothed

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    fname = "sdmf_pyxelmask.h5"
    fid = h5py.File(fname, "r")
    ds_orbits = fid["orbits"]
    idx = np.argsort(ds_orbits[:])
    orbits = ds_orbits[:][idx]
    print(idx)
    ds_combi = fid["combined"]
    num_orbits = idx.size
    a = np.empty((num_orbits,1024), dtype=np.float)
    for i_orbit in range(num_orbits):
        id_ = idx[i_orbit]
        a[i_orbit,:] = ds_combi[id_,:]
    fid.close()

    # print(a.shape)
    # plt.cla()
    # plt.imshow(a)
    # plt.show()

    # fname_out = "orbital_dbqm.h5"
    # fid_out = h5py.File(fname_out, "w")
    # ds = fid_out.create_dataset("data", a.shape, dtype=np.float)
    # ds[:,:] = a
    # fid_out.close()
    # np.savetxt("orbital_dbqm.csv", a, delimiter=",")
    # scipy.misc.imsave('orbital_dbqm.png', a*255)

    s = smooth(orbits, a)

    print(s.shape)
    pixnr = 450
    orig = a[:,pixnr]
    new = s[:,pixnr]
    ax = np.array(range(num_orbits))
    print(orig)
    print(new)
    plt.cla()
    # plt.imshow(s)
    plt.plot(ax, orig, 'bo', ax, new, 'ro')
    plt.show()
