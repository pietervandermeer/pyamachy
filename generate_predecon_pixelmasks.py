#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from pixelquality import PixelQuality
from sciamachy import decon_start_orbits
#from envisat import last_orbit

if __name__ == "__main__":
    pq = PixelQuality()

    # orbits = decon_start_orbits
    # orbits.append(52000) # add end of mission to have a full list of intervals "between" decontaminations
    #orbits = range(6200,52000)
    orbits = range(5700,8300)

    for orbit in orbits:
#        if orbit < 6000:
#            continue # ehehehe, let's just skip that one for now..
#        orbit -= 200 # just before start of decontamination
        print(orbit)
        pq.calculate(orbit)
        pq.write(directory=".")
        pq.write_ascii(directory="qualities", all_figures=False)
