#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
average pixel quality over entire mission
"""

from __future__ import print_function, division

import h5py
import numpy as np
from pixelquality import PixelQuality
from sciamachy import decon_start_orbits

if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    fname = "sdmf_pyxelmask.h5"
    fid = h5py.File(fname, "r")
    ds_orbit = fid["orbits"]
    ds_combined = fid["combined"]
    ds_invalid = fid["invalid"]
#    ds_combined = fid["combined"]
#    ds_combined = fid["combined"]

    average_quality = np.zeros(1024, dtype=np.float64)
    average_invalid = np.zeros(1024, dtype=np.float64)

    #np.nan_to_num()
    #np.mean(ds_combined[:,:], axis=0) # nan problems

    print(ds_combined.shape)

    i = 0
    n = 0
    for orbit in ds_orbit[:]:
        #print(ds_combined[i,:])
        num_nan = np.sum(np.isnan(ds_combined[i,:]))
        if num_nan == 0:
            #print(orbit)
            average_quality += ds_combined[i,:]
            average_invalid += ds_invalid[i,:]
            n += 1
        else:
            print(orbit)
        i += 1

    average_quality /= n
    average_invalid /= n
    print("rows summed =", n)

    #print(average_quality)
    for elem in average_quality:
        print(elem)
    # for elem in average_invalid:
    #     print(elem)
