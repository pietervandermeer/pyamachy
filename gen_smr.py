#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Generates smr databases in various levels of calibration: none (raw), non-linearity&memory (12), nlin&mem+dark (123).
"""

import argparse

import generate_smr
from envisat import parseOrbitList

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
        orbitrange = "4152-53200"

    generate_smr.main("db -cnone --orbit="+orbitrange+" sdmf_smr_raw.h5")
    generate_smr.main("db -c1,2 --orbit="+orbitrange+" sdmf_smr_12.h5")
    generate_smr.main("db -c1,2,3 --orbit="+orbitrange+" sdmf_smr_123.h5")

