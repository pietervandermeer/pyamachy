#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Smooths channel 8 pixel quality mask.
"""

from __future__ import print_function, division

import numpy as np 
import h5py
import matplotlib.pyplot as plt

from pixelquality_module import create_figure_dset, create_mask_dset
from mask import Mask

#-- globals --------------------------------------------------------------------

input_mask_fname = "sdmf_pyxelmask.h5"

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
    #print(num_orbits)

    #
    # sort data rows on orbit number
    #

    idx = np.argsort(orbits)
    #print("sorted orbits", orbits[idx])
    #print("raw data", data)
    data = data[idx,:]
    #print("indexed data", data)

    for i_orbit in range(num_orbits):
        upper = i_orbit + winsize/2
        lower = i_orbit - winsize/2
        if upper > num_orbits:
            upper = num_orbits
        if lower < 0:
            lower = 0
        window = np.arange(lower, upper, dtype=np.int)
        #print(window)
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
    """
    Load sdmf (orbital) pixelmask figure from database.

    Parameters
    ----------

    start_orbit : int
        start orbit
    end_orbit : int
        end orbit
    fig_name : str
        figure name

    Returns
    -------

    orbits32 : array, dtype='u2'
        orbit list 
    a : array, dtype=float
        figure data
    """
    fname = input_mask_fname
    fid = h5py.File(fname, "r")
    ds_orbits = fid["orbits"]
    orbits = ds_orbits[:]

    idx = np.argsort(orbits[(orbits >= start_orbit) & (orbits <= end_orbit)])
    #print(idx)

    orbits32 = orbits[(orbits >= start_orbit) & (orbits <= end_orbit)][idx]
    #print("orbits32", orbits32)

    orbit_range = end_orbit - start_orbit + 1
    ds_combi = fid[fig_name]
    num_orbits = idx.size
    a = np.empty((num_orbits,1024), dtype=np.float)
    i_row = 0
    for orbit in orbits32:
        i_orbit = np.where(orbit == ds_orbits[:])[0][0]
        #print(ds_combi.dtype)
        if ds_combi.dtype != np.bool:
            a[i_row,:] = np.nan_to_num(ds_combi[i_orbit,:])
        else:
            a[i_row,:] = ds_combi[i_orbit,:]
        i_row += 1

        #print("i_orbit", i_orbit, "row", a[i_row,:])
    fid.close()
    return orbits32, a 

def compare_smoothmasks():
    fname = input_mask_fname #"sdmf_pyxelmask.h5"
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
    import argparse
    import warnings
    warnings.simplefilter("error") # warnings to errors
    from envisat import parseOrbitList
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    num_pixels = 1024

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', dest='input_fname', type=str, help="input (orbital) mask file name")
    parser.add_argument('-o', '--output', dest='output_fname', type=str, 
                        help="output file name")
    parser.add_argument('-c', action='store_true', dest="sdmf30_compat")
    parser.add_argument('--config', dest='config_file', type=file, 
                        default='default3.2.cfg')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', 
                        version='%(prog)s 0.1')
    parser.add_argument('-P', '--path', dest='path')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--orbitrange', default='all', 
                        help='sets orbit range f.e. "43000-44000", "all"', type=parseOrbitList)
    parser.add_argument('-p', '--pixnr', action='store', type=int, default=None,
                        dest='pixnr', help="pixel number to be examined [0..1023]")
    args = parser.parse_args()

    #
    # set defaults
    #

    start_orbit = 4152
    end_orbit = 53000
    output_fname = "sdmf_smooth_pyxelmask.h5"

    #
    # handle command line arguments
    #

    sdmf30_compat = args.sdmf30_compat
    verbose = args.verbose
    orbitrange = args.orbitrange
    if orbitrange is not None:
        start_orbit = orbitrange[0]
        end_orbit = orbitrange[1]
    if args.input_fname is not None:
        input_mask_fname = args.input_fname
    if args.output_fname is not None:
        output_fname = args.output_fname

    print("output_fname =", output_fname)
    print("input_mask_fname =", input_mask_fname)
    print("sdmf30_compat =", sdmf30_compat)
    print("orbitrange =", start_orbit, end_orbit)

    #
    # load all pixel quality data, figure by figure
    #

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
    print("orbits.size", orbits.size, "inv_data.shape", inv_data.shape)

    #
    # smooth the figures
    #

    print("smoothing..")
    inv_smooth = 1.0 - smooth(orbits, inv_data) # input was boolean flags (1: bad, 0: good).. so inverted
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

    if sdmf30_compat:
        combined = inv_smooth * chi_smooth * (noise_smooth*noise_smooth) * sat_smooth * sun_smooth * wls_smooth
    else:
        # TODO: use weights from config file
        combined = inv_smooth * (0.5*darkerr_smooth + 0.5*darkres_smooth) * sat_smooth * sun_smooth * wls_smooth

    combined_flag = np.nan_to_num(combined) < 0.1

    print("combined.")

    #
    # write smooth and combined data
    #

    print("writing hdf5..")

    fid = h5py.File(output_fname, "w")
    fid.attrs['sdmf30_compat']=sdmf30_compat

    ds = fid.create_dataset("orbits", orbits.shape, dtype=np.int, data=orbits)
    ds.attrs["long_name"] = np.string_("Absolute orbit number")

    ds = create_figure_dset(fid, "combined", dims=combined.shape, data=combined)
    ds.attrs["long_name"] = np.string_("Smooth combined figure (channel 8)")
    ds.attrs["units"] = np.string_("-")
    ds.attrs["description"] = np.string_("Figure [0.0..1.0] that combines all the other figures to estimate a total quality for the pixel (channel 8).")

    ds = create_mask_dset(fid, "combinedFlag", dims=combined_flag.shape, data=combined_flag)
    ds.attrs["long_name"] = np.string_("Smooth combined flag (channel 8)")
    ds.attrs["units"] = np.string_("-")
    ds.attrs["description"] = np.string_("""Boolean flag that indicates a pixel as good (0), or bad (1). 
                                         It combines all criteria (= `combined' figure thesholded).""")

    ds = create_figure_dset(fid, "invalid", dims=inv_smooth.shape, data=inv_smooth)
    ds.attrs["long_name"] = np.string_("Smooth invalid flag (channel 8)")
    ds.attrs["units"] = np.string_("-")
    ds.attrs["description"] = np.string_("Boolean flag that indicates if the dark fit has succeeded (1) or failed (0).")

    ds = create_figure_dset(fid, "saturation", dims=sat_smooth.shape, data=sat_smooth)
    ds.attrs["long_name"] = np.string_("Smooth saturation quality figure (channel 8)")
    ds.attrs["units"] = np.string_("-")
    ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if the dark current is fully saturated (1.0), 
                                         not saturated (0.0), or somewhere inbetween.
                                         This figure is only nonzero if darks are closer than a few thousand BU removed from saturation.""")

    ds = create_figure_dset(fid, "sunResponse", dims=sun_smooth.shape, data=sun_smooth)
    ds.attrs["long_name"] = np.string_("Smooth sun response quality figure (channel 8)")
    ds.attrs["units"] = np.string_("-")
    ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that gives relative deviation from expected Sun-over-ESM diffuser (state 62) response.
                                         For instance, a 10%% response is represented as 0.1, a 1000%% response is too.
                                         """)

    ds = create_figure_dset(fid, "wlsResponse", dims=wls_smooth.shape, data=wls_smooth)
    ds.attrs["long_name"] = np.string_("Smooth WLS Response quality figure (channel 8)")
    ds.attrs["units"] = np.string_("-")
    ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that gives relative deviation from expected WLS (White Light Source) response.
                                         For instance, a 10%% response is represented as 0.1, a 1000%% response is too.
                                         """)
    if sdmf30_compat:
        ds = create_figure_dset(fid, "chisquare3.0", dims=chi_smooth.shape, data=chi_smooth)
        ds.attrs["long_name"] = np.string_("Smooth dark Chi^2 from SDMF3.0 (channel 8)")
        ds.attrs["units"] = np.string_("-")
        ds.attrs["description"] = np.string_("""Chi-squared from the dark signal fit for channel 8. 
                                             This is SDMF3.0 data, which may be used for backward compatibility.""")

        ds = create_figure_dset(fid, "noise", dims=noise_smooth.shape, data=noise_smooth)
        ds.attrs["long_name"] = np.string_("Smooth noise criterion from SDMF3.0 (channel 8)")
        ds.attrs["units"] = np.string_("-")
        ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if a pixel exceeds its expected noise 
                                             (from beginning of mission) by a large margin. 
                                             This is SDMF3.0 data, which may be used for backward compatibility.""")
    else:
        ds = create_figure_dset(fid, "darkError", dims=darkerr_smooth.shape, data=darkerr_smooth)
        ds.attrs["long_name"] = np.string_("Smooth dark correction error quality figure (channel 8)")
        ds.attrs["units"] = np.string_("-")
        ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if the dark correction error (darks corrected by vardark) 
                                             is too big (1.0), none (0.0), or somewhere inbetween.""")

        ds = create_figure_dset(fid, "darkResidual", dims=darkres_smooth.shape, data=darkres_smooth)
        ds.attrs["long_name"] = np.string_("Smooth dark fit residual quality figure (channel 8)")
        ds.attrs["units"] = np.string_("-")
        ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if the fit residual is too big (1.0), none (0.0), or somewhere inbetween.""")

    fid.close()

    print("written.")
