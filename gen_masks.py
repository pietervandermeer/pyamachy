#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Generate all pixelmasks: orbital, and smoothed, SDMF3.0 compatible and SDMF3.2.
"""

import argparse

from envisat import parseOrbitList
import generate_pyxelmask
import generate_smooth_pyxelmask

if __name__ == "__main__":

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', dest='loglevel', type=str, 
                        choices=("DEBUG", "INFO", "WARN", "ERROR", "FATAL"), 
                        default="INFO", help="logging level")
    parser.add_argument('-V', '--version', action='version', 
                        version='%(prog)s 0.1')
    parser.add_argument('-P', '--path', dest='path', default=".")
    parser.add_argument('--orbitrange', default='all', 
                        help='sets orbit range f.e. "43000-44000", "all"', type=parseOrbitList)
    args = parser.parse_args()

    if args.orbitrange is None:
        orbitrange = 4152,53200

    general_args = " --log "+args.loglevel+" --orbitrange="+orbitrange

    orbitalmask30_fname = path+"/sdmf_pyxelmask30.h5"
    orbitalmask32_fname = path+"/sdmf_pyxelmask32.h5"

    generate_pyxelmask.main("-c -o "+orbitalmask30_fname+general_args)
    generate_pyxelmask.main("-o "+orbitalmask32_fname+general_args)

    smoothmask30_fname = path+"/sdmf_smooth_pyxelmask30.h5"
    smoothmask32_fname = path+"/sdmf_smooth_pyxelmask32.h5"

    generate_smooth_pyxelmask.main("-c -i "+orbitalmask30_fname+" -o "+smoothmask30_fname+general_args)
    generate_smooth_pyxelmask.main(" -i "+orbitalmask32_fname+" -o "+smoothmask32_fname+general_args)
