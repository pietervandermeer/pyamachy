#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Interpolates vardark fit parameters between monthly orbits.
This uses a monthly_fits*.h5 database as input.
Interpolated parameters are generated and stored for every orbit.

Note that it will also extrapolate beyond the first and last orbit when necessary.
Take care that a partial monthly_fits database will be extrapolated over the entire mission.
"""

from __future__ import print_function, division

import warnings
warnings.simplefilter("error") # warnings to errors
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy import array
import h5py

from envisat import parseOrbitList
from sciamachy import orbitfilter
from vardark import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
from fit_monthlies import fit_monthlies

#-- globals --------------------------------------------------------------------

n_pix = 1024
first_orbit = 5000
last_orbit = 6000

#-- functions ------------------------------------------------------------------

def extrap1d(interpolator):
    """
    Extrapolate outside of the original data range of the specified interpolator
    """

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

def interpolate_monthlies(db_out_name, db_in_name):
    monthlies_dbname = db_in_name
    if not h5py.is_hdf5(monthlies_dbname):
        print(monthlies_dbname+" not present. generating...")
        fit_monthlies(monthlies_dbname)
    fin = h5py.File(monthlies_dbname, "r")
    monthlies = fin["orbits"]
    #n_orbits = monthlies[-1] - monthlies[0]
    #orblist = np.arange(monthlies[0], monthlies[-1])
    n_orbits = 51996 - monthlies[0] # also include last orbits, beyond last monthly calibration orbit
    print("n_orbits=",n_orbits)
    orblist = np.arange(monthlies[0], 51996)
    aos_dset = fin["aos"]
    off_dset = fin["off"]
    amps_dset = fin["amps"]
    phase_dset = fin["channel_phase"]
    amp2_dset = fin["amp2"]
    err_aos_dset = fin["err_aos"]
    err_off_dset = fin["err_off"]
    err_amps_dset = fin["err_amps"]
    err_phase_dset = fin["err_channel_phase"]
    err_amp2_dset = fin["err_amp2"]

    # create output hdf5

    fout = h5py.File(db_out_name, "w")
    iorbits = fout.create_dataset("orbits", (n_orbits,), dtype='i')
    iorbits.attrs["long_name"] = np.string_("Absolute orbit numbers.")
    iorbits.attrs["units"] = np.string_("-")
    iorbits.attrs["description"] = np.string_("Absolute orbit numbers.")

    iphases = fout.create_dataset("phases", (n_orbits,2), dtype='f')
    iphases.attrs["long_name"] = np.string_("Phase shifts of channel 8 orbital variation")
    iphases.attrs["units"] = np.string_("6.02e3 s")
    iphases.attrs["description"] = np.string_("Orbit phase shifts wrt eclipse of fundamental and first harmonic [0..1] per orbit.")

    iaos = fout.create_dataset("aos", (n_orbits,n_pix), dtype='f')
    iaos.attrs["long_name"] = np.string_("analog offset (channel 8)")
    iaos.attrs["units"] = np.string_("BU")
    iaos.attrs["description"] = np.string_("SCIA analog offset of channel 8")

    iamps = fout.create_dataset("amps", (n_orbits,n_pix), dtype='f')
    iamps.attrs["long_name"] = np.string_("Fundamental amplitude, channel 8")
    iamps.attrs["units"] = np.string_("BU/s")
    iamps.attrs["description"] = np.string_("Amplitude of fundamental sinusoid, channel 8, per orbit per pixel.")

    iamp2 = fout.create_dataset("amp2", (n_orbits,), dtype='f')
    iamp2.attrs["long_name"] = np.string_("Amplitude of first harmonic, channel 8")
    iamp2.attrs["units"] = np.string_("")
    iamp2.attrs["description"] = np.string_("Amplitude of first harmonic, relative to fundamental amplitude [0..1], channel average.")

    i_err_phases = fout.create_dataset("err_phases", (n_orbits,n_pix,2), dtype='f')
    i_err_phases.attrs["long_name"] = np.string_("Error in phase shifts of channel 8 orbital variation")
    i_err_phases.attrs["units"] = np.string_("6.02e3 s")
    i_err_phases.attrs["description"] = np.string_("Error in orbit phase shifts wrt eclipse of fundamental and first harmonic [0..1] per orbit per pixel.")

    i_err_aos = fout.create_dataset("err_aos", (n_orbits,n_pix), dtype='f')
    i_err_aos.attrs["long_name"] = np.string_("Error in analog offset (channel 8)")
    i_err_aos.attrs["units"] = np.string_("BU")
    i_err_aos.attrs["description"] = np.string_("Error in SCIA analog offset of channel 8, per orbit per pixel")

    i_err_amps = fout.create_dataset("err_amps", (n_orbits,n_pix), dtype='f')
    i_err_amps.attrs["long_name"] = np.string_("Error in fundamental amplitude, channel 8.")
    i_err_amps.attrs["units"] = np.string_("BU/s")
    i_err_amps.attrs["description"] = np.string_("Error in amplitude of fundamental sinusoid, channel 8, per orbit per pixel.")

    i_err_amp2 = fout.create_dataset("err_amp2", (n_orbits,n_pix), dtype='f')
    i_err_amp2.attrs["long_name"] = np.string_("Error in amplitude of first harmonic, channel 8")
    i_err_amp2.attrs["units"] = np.string_("")
    i_err_amp2.attrs["description"] = np.string_("Error in amplitude of first harmonic, relative to fundamental amplitude [0..1], per orbit per pixel.")

    fun_phase1 = extrap1d(interp1d(monthlies, phase_dset[:,0]))
    fun_phase2 = extrap1d(interp1d(monthlies, phase_dset[:,1]))
    fun_amp2 = extrap1d(interp1d(monthlies, amp2_dset[:]))
    fun_aos = extrap1d(interp1d(monthlies, aos_dset[:,:], axis=0))
    fun_amps = extrap1d(interp1d(monthlies, amps_dset[:,:], axis=0))

    fun_err_phase1 = extrap1d(interp1d(monthlies, err_phase_dset[:,:,0], axis=0))
    fun_err_phase2 = extrap1d(interp1d(monthlies, err_phase_dset[:,:,1], axis=0))
    fun_err_amp2 = extrap1d(interp1d(monthlies, err_amp2_dset[:,:], axis=0))
    fun_err_aos = extrap1d(interp1d(monthlies, err_aos_dset[:,:], axis=0))
    fun_err_amps = extrap1d(interp1d(monthlies, err_amps_dset[:,:], axis=0))

    iorbits[:] = orblist

    iaos[:,:] = fun_aos(orblist)
    iamps[:,:] = fun_amps(orblist)
    iamp2[:] = fun_amp2(orblist)
    iphases[:,0] = fun_phase1(orblist)
    iphases[:,1] = fun_phase2(orblist)

    i_err_aos[:,:] = fun_err_aos(orblist)
    i_err_amps[:,:] = fun_err_amps(orblist)
    i_err_amp2[:,:] = fun_err_amp2(orblist)
    i_err_phases[:,:,0] = fun_err_phase1(orblist)
    i_err_phases[:,:,1] = fun_err_phase2(orblist)

    # write out and close
    fout.close()
    # close input
    fin.close()

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser, ArgumentTypeError

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', dest='input_fname', type=str, default='monthly_fits_long.h5')
    parser.add_argument('-o', '--output', dest='output_fname', type=str, default='interpolated_monthlies_long.h5')
    parser.add_argument('--config', dest='config_file', type=file, 
                        default='default3.1.cfg')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', 
                        version='%(prog)s 0.1')
    parser.add_argument('--orbitrange', default='all', 
                        help='sets orbit range f.e. "43000-44000", "all"', type=parseOrbitList)
    # parameters used especially for the plot
    parser.add_argument('-s', '--short', action='store_true', 
                        dest='shortMode', help="compute product for short PET's instead of long ones.")
    args = parser.parse_args()

    #
    # create default output name
    #

    if args.shortMode:
        input_fname = 'monthly_fits_short.h5'
        output_fname = "interpolated_monthlies_short.h5"
    else:
        input_fname = 'monthly_fits_long.h5'
        output_fname = "interpolated_monthlies_long.h5"

    print(input_fname, "-> interpolate ->", output_fname)

    interpolate_monthlies(output_fname, input_fname)
