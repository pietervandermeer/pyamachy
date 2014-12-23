#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Fits vardark parameters to monthly calibration orbits and writes to database monthly_fits_*.h5. 
NOTE: file will be truncated!
"""

from __future__ import print_function, division

import numpy as np
import h5py
import re
from scipy.interpolate import interp1d

from envisat import parseOrbitList
from sciamachy import orbitfilter, get_darkstateid, petcorr, n_chanpix
from vardark import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m

n_pix = n_chanpix
pixnr = 597
first_orbit = 1990
last_orbit = 55000

#-- functions ------------------------------------------------------------------

def is_blacklisted(orbit):
    """
    Unfittable monthly calibration orbits -> ditch 'em!
    """
    return (orbit == 2216) or (orbit == 3921)

def no_monthly(orbit, old_monthly, last_orbit):
    """
    Returns if specified orbit is a new monthly calibration orbit and in range.
    Some monthlies are corrupted and hence are blacklisted.
    """
    return (orbit == old_monthly) or (orbit > last_orbit)

def fit_monthlies(db_out_name, short=False):
    ofilt = orbitfilter()
    if short:
        ad = AllDarks([1.0, 0.5, 0.125, 0.0625]) 
    else:
        ad = AllDarks([1.0, 0.5, 0.125])

    # compute amount of monthlies
    n_monthlies = 0
    old_monthly = ofilt.get_next_monthly(first_orbit)
    while True:
        monthly = ofilt.get_next_monthly(old_monthly)
        if no_monthly(monthly, old_monthly, last_orbit):
            break
        old_monthly = monthly
        if not is_blacklisted(monthly):
            n_monthlies += 1

    f = h5py.File(db_out_name, "w")
    orblist = f.create_dataset("orbits", (n_monthlies,), dtype='i')
    aos_dset = f.create_dataset("aos", (n_monthlies,n_pix), dtype='f')
    off_dset = f.create_dataset("off", (n_monthlies,n_pix), dtype='f')
    amps_dset = f.create_dataset("amps", (n_monthlies,n_pix), dtype='f')
    phase_dset = f.create_dataset("channel_phase", (n_monthlies,2), dtype = 'f')
    amp2_dset = f.create_dataset("amp2", (n_monthlies,), dtype='f')
    err_aos_dset = f.create_dataset("err_aos", (n_monthlies,n_pix), dtype='f')
    err_off_dset = f.create_dataset("err_off", (n_monthlies,n_pix), dtype='f')
    err_amps_dset = f.create_dataset("err_amps", (n_monthlies,n_pix), dtype='f')
    err_phase_dset = f.create_dataset("err_channel_phase", (n_monthlies,n_pix,2), dtype = 'f')
    err_amp2_dset = f.create_dataset("err_amp2", (n_monthlies,n_pix), dtype='f')

    old_monthly = ofilt.get_next_monthly(first_orbit)
    i_monthly = 0
    while True:
        monthly = ofilt.get_next_monthly(old_monthly)
        if no_monthly(monthly, old_monthly, last_orbit):
            break
        if not is_blacklisted(monthly):
            # compute monthly fit, to obtain analog offset
            print("compute monthly fit:", monthly, "..")
            channel_phase1, channel_phase2, aos, off, amps, channel_amp2, trends, errors = fit_monthly(ad, monthly, short=short, give_errors=True)
            aos_dset[i_monthly,:] = aos
            off_dset[i_monthly,:] = off
            amps_dset[i_monthly,:] = amps
            phase_dset[i_monthly,:] = np.array([channel_phase1, channel_phase2])
            amp2_dset[i_monthly] = channel_amp2
            err_aos_dset[i_monthly,:] = errors["aos"]
            err_off_dset[i_monthly,:] = errors["off"]
            err_amps_dset[i_monthly,:] = errors["amps"]
            err_phase_dset[i_monthly,:,0] = errors["phase1"]
            err_phase_dset[i_monthly,:,1] = errors["phase2"]
            err_amp2_dset[i_monthly,:] = errors["amps2"]
            orblist[i_monthly] = monthly
            i_monthly += 1
        old_monthly = monthly

    print("done.")

    f.close()

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser, ArgumentTypeError

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output', dest='output_fname', type=str)
    parser.add_argument('--config', dest='config_file', type=file, 
                        default='default3.1.cfg')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', 
                        version='%(prog)s 0.1')
    parser.add_argument('--orbitrange', help='sets orbit range f.e. "43000-44000", "all"', type=parseOrbitList)
    # parameters used especially for the plot
    parser.add_argument('-s', '--short', action='store_true', 
                        dest='shortMode', help="compute product for short PET's instead of long ones.")
    args = parser.parse_args()

    #
    # create default output name
    #

    if args.shortMode:
        output_fname = "monthly_fits_short.h5"
    else:
        output_fname = "monthly_fits_long.h5"

    if args.output_fname is not None:
        output_fname = args.output_fname

    if args.orbitrange is not None:
        first_orbit = args.orbitrange[0]
        last_orbit = args.orbitrange[1]

    fit_monthlies(output_fname, short=args.shortMode)
