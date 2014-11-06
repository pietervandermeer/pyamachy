from __future__ import print_function, division

import numpy as np
import h5py
import re
from scipy.interpolate import interp1d

from sciamachy_module import orbitfilter, get_darkstateid, petcorr, n_chanpix
from vardark_module import AllDarks, trending_phase, fit_monthly, fit_eclipse_orbit
from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m

n_pix = n_chanpix
pixnr = 597
first_orbit = 5000
last_orbit = 55000

#-- functions ------------------------------------------------------------------

def parseOrbitList(str):
    msg1 = "'" + str + "' is not a range or number." \
        + " Expected forms like '20000-25000' or '20000'."
    msg2 = "'" + str + "' is not valid orbit number."

    if str.lower() == 'all':
        return None

    m = re.match(r'(\d+)(?:-(\d+))?$', str)
    if not m:
        raise ArgumentTypeError( msg1 )
    v1 = int(m.group(1))
    if m.group(2):
        v2 = int(m.group(2))
        if v1 < 1 or v2 > 100000:
            raise ArgumentTypeError( msg2 )
        return (v1, v2)
    else:
        return v1

def fit_monthlies(db_out_name, short=False):
    ofilt = orbitfilter()
    if short:
        ad = AllDarks([0.125, 0.06250])
    else:
        ad = AllDarks([1.0, 0.5])

    # compute amount of monthlies
    n_monthlies = 0
    old_monthly = ofilt.get_next_monthly(first_orbit)
    while True:
        monthly = ofilt.get_next_monthly(old_monthly)
        if (monthly == old_monthly) or (monthly > last_orbit):
            break
        old_monthly = monthly
        n_monthlies += 1

    f = h5py.File(db_out_name, "w")
    orblist = f.create_dataset("orbits", (n_monthlies,), dtype='i')
    aos_dset = f.create_dataset("aos", (n_monthlies,n_pix), dtype='f')
    lcs_dset = f.create_dataset("lcs", (n_monthlies,n_pix), dtype='f')
    amps_dset = f.create_dataset("amps", (n_monthlies,n_pix), dtype='f')
    phase_dset = f.create_dataset("channel_phase", (n_monthlies,2), dtype = 'f')
    amp2_dset = f.create_dataset("amp2", (n_monthlies,), dtype='f')

    old_monthly = ofilt.get_next_monthly(first_orbit)
    i_monthly = 0
    while True:
        monthly = ofilt.get_next_monthly(old_monthly)
        if (monthly == old_monthly) or (monthly > last_orbit):
            break
        # compute monthly fit, to obtain analog offset
        print("compute monthly fit:", monthly, "..")
        channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends = fit_monthly(ad, monthly)
        aos_dset[i_monthly,:] = aos
        lcs_dset[i_monthly,:] = lcs
        amps_dset[i_monthly,:] = amps
        phase_dset[i_monthly,:] = np.array([channel_phase1, channel_phase2])
        amp2_dset[i_monthly] = channel_amp2
        old_monthly = monthly
        orblist[i_monthly] = monthly
        i_monthly += 1

    print("done.")

    f.close()

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser, ArgumentTypeError

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=
             'fits vardark parameters to monthly calibration orbits and writes to database monthly_fits_*.h5. NOTE: file will be truncated!')
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
