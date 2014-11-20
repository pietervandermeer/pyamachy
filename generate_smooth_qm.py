from __future__ import print_function, division

import numpy as np 
import h5py
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import scipy.misc
import warnings
warnings.simplefilter("error") # warnings to errors

from mask import Mask

#-- globals --------------------------------------------------------------------


#-- functions ------------------------------------------------------------------

def smooth(orbits, data, winsize=50):
    """
    Smooths multiple pixels over time, using specified window. 
    Missing orbits are not taken into account: 
    f.e. orbit 1000 can be orbit 2000's neighbor when the range inbetween is missing.

    Parameters
    ----------

    orbits : 1d numpy array, int
        array of orbits, every entry corresponds to one row in `data' (not used yet)
    data : 2d array, shape (orbits,pixels), dtype flexible
        any type of data: bool, int or float 
    winsize : int, optional
        smoothing window size in number of orbits

    Returns
    -------

    smoothed version of `data'. same dimensions, dtype float. 
    """

    smoothed = np.empty(data.shape, dtype=np.float64)
    num_orbits = data.shape[0]

    for i_orbit in range(num_orbits):
        if i_orbit < winsize/2:
            window = np.arange(0, i_orbit+winsize/2, dtype=np.int)
        elif i_orbit >= num_orbits - winsize/2:
            window = np.arange(i_orbit-winsize/2, num_orbits, dtype=np.int)
        else:
            window = np.arange(i_orbit-winsize/2, i_orbit+winsize/2, dtype=np.int)
        smoothed[i_orbit, :] = np.mean(data[window,:], axis=0)

    return smoothed

def plot(orbits_f, flt, orbits_b, bin, pixnr):
    plt.cla()
#    fig = plt.figure()
    fig = plt.subplot(111)
    fig.set_ylim([-.1,1.1])
    fig.set_title("pixel "+str(pixnr))

    print(np.min(orbits_f), np.max(orbits_f), np.min(orbits_b), np.max(orbits_b))
    s_bin = smooth(orbits_b, bin)
    bin_slice = (1.0 - bin[:,pixnr].flatten())
    sbin_slice = (1.0 - s_bin[:,pixnr].flatten())
    ax = np.arange(bin_slice.size)
#    plt.plot(ax, bin_slice, 'ro', ax, sbin_slice, 'r-')
    plt.plot(orbits_b, bin_slice, 'r-')

    s_flt = smooth(orbits_f, flt)
    flt_slice = flt[:,pixnr].flatten()
    sflt_slice = s_flt[:,pixnr].flatten()
    ma = np.max(flt_slice)
    if ma > 0:
        flt_slice /= ma
        sflt_slice /= ma
    ax = np.arange(flt_slice.size)
    plt.plot(orbits_f, flt_slice, 'bo', orbits_f, sflt_slice, 'g-')

    plt.show()
    return

def load_sdmf32_figure(start_orbit, end_orbit, fig_name):
    fname = "sdmf_pyxelmask.h5"
    fid = h5py.File(fname, "r")
    ds_orbits = fid["orbits"]
    idx = np.argsort(ds_orbits[:])
    orbits32 = ds_orbits[:][idx]
    start_orbit = np.min(orbits32)
    stop_orbit = np.max(orbits32)
    orbit_range = stop_orbit - start_orbit + 1
    ds_combi = fid[fig_name]
    num_orbits = idx.size
    a = np.empty((num_orbits,1024), dtype=np.float)
    for i_orbit in range(num_orbits):
        id_ = idx[i_orbit]
        a[i_orbit,:] = ds_combi[id_,:]
    fid.close()
    return orbits32, a 

def compare_smoothmasks():
    fname = "sdmf_pyxelmask.h5"
    fid = h5py.File(fname, "r")
    ds_orbits = fid["orbits"]
    idx = np.argsort(ds_orbits[:])
    orbits32 = ds_orbits[:][idx]
    start_orbit = np.min(orbits32)
    stop_orbit = np.max(orbits32)
    orbit_range = stop_orbit - start_orbit + 1
    #print(idx)
    ds_combi = fid["combined"]
    num_orbits = idx.size
    a = np.empty((num_orbits,1024), dtype=np.float)
    for i_orbit in range(num_orbits):
        id_ = idx[i_orbit]
        orbit = orbits32[i_orbit]
#        a[orbit-start_orbit,:] = ds_combi[id_,:]
        a[i_orbit,:] = ds_combi[id_,:]
    fid.close()
    print("3.2 read")

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

    m = Mask()
    a_binary = np.ones((1400,1024), dtype=np.bool)
    i_orbit = 0
    orbits30 = np.empty(1400, dtype=np.int)
    for orbit in range(42000,43400):
        #print(orbit)
        try:
            m.load_sdmf30_crit(orbit, "combined", smooth=True)
        except:
            continue
        orbits30[i_orbit] = orbit
        a_binary[i_orbit,:] = m.mask
        i_orbit += 1
    print("3.0 read")

    orbits30 = orbits30[0:i_orbit]
    a_binary = a_binary[0:i_orbit,:]

    for pixnr in range(1024):
        plot(orbits32, a, orbits30, a_binary, pixnr)

    return

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    # TODO: commandline interface
    sdmf30_compat = True
    num_pixels = 1024

    #
    # load all pixel quality data, figure by figure
    #

    start_orbit = 5593
    end_orbit = 53000

    orbits, inv_data = load_sdmf32_figure(start_orbit, end_orbit, "invalid")
    print("read invalid")
    orbits, sat_data = load_sdmf32_figure(start_orbit, end_orbit, "saturation")
    print("read saturation")
    orbits, sun_data = load_sdmf32_figure(start_orbit, end_orbit, "sunResponse")
    print("read sunResponse")
    orbits, wls_data = load_sdmf32_figure(start_orbit, end_orbit, "wlsResponse")
    print("read wlsResponse")

    if sdmf30_compat:
        orbits, chi_data = load_sdmf32_figure(start_orbit, end_orbit, "chisquare3.0")
        print("read chisquare3.0 (SDMF3.0 compatibility)")
        orbits, noise_data = load_sdmf32_figure(start_orbit, end_orbit, "noise")
        print("read noise (SDMF3.0 compatibility)")
    else:
        orbits, darkerr_data = load_sdmf32_figure(start_orbit, end_orbit, "darkError")
        print("read darkError")
        orbits, darkres_data = load_sdmf32_figure(start_orbit, end_orbit, "darkResidual")
        print("read darkResidual")

    num_orbits = orbits.size

    #
    # smooth the figures
    #

    print("smoothing..")
    inv_smooth = smooth(orbits, sat_data)
    sat_smooth = smooth(orbits, sat_data)
    sun_smooth = smooth(orbits, sun_data, winsize=100)
    wls_smooth = smooth(orbits, wls_data, winsize=500)
    if sdmf30_compat:
        chi_smooth = smooth(orbits, chi_data)
        noise_smooth = smooth(orbits, noise_data)
    else:
        darkerr_smooth = smooth(orbits, darkerr_data)
        darkres_smooth = smooth(orbits, darkres_data)
    print("smoothed.")

    #
    # combine smoothed figures
    #

    print("combining smoothed figures..")
 
    combined = np.empty(sat_data.shape, dtype=np.float64)

    for i_orbit in range(num_orbits):
        if sdmf30_compat:
            combined = inv_smooth * chi_smooth * (noise_smooth**2) * sat_smooth * sun_smooth * wls_smooth
        else:
            # TODO: use weights from config file
            combined = inv_smooth * (0.5*darkerr_smooth + 0.5*darkres_smooth) * sat_smooth * sun_smooth * wls_smooth

    print("combined.")

    #
    # write smooth and combined data
    #

    print("writing hdf5..")

    fid = h5py.File("sdmf_smooth_pyxelmask.h5", "w")

    fid.create_dataset("orbits", orbits.shape, dtype=np.int, data=orbits)
    fid.create_dataset("invalid", inv_smooth.shape, dtype=np.float, data=inv_smooth)
    fid.create_dataset("saturation", sat_smooth.shape, dtype=np.float, data=sat_smooth)
    fid.create_dataset("sunResponse", sun_smooth.shape, dtype=np.float, data=sun_smooth)
    fid.create_dataset("wlsResponse", wls_smooth.shape, dtype=np.float, data=wls_smooth)
    if sdmf30_compat:
        fid.create_dataset("chisquare3.0", chi_smooth.shape, dtype=np.float, data=chi_smooth)
        fid.create_dataset("noise", noise_smooth.shape, dtype=np.float, data=noise_smooth)
    else:
        fid.create_dataset("darkError", darkerr_smooth.shape, dtype=np.float, data=darkerr_smooth)
        fid.create_dataset("darkResidual", darkres_smooth.shape, dtype=np.float, data=darkres_smooth)

    fid.close()

    print("written.")