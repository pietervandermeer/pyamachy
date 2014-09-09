# -*- coding: iso-8859-1 -*-
"""
Contains functionality related to ENVISAT satellite operations.
"""

import numpy
import time
from datetime import date

# returns MJD dates (starting 1-1-1) for given ENVISAT absolute orbit numbers
def convert_orbit_to_jd(orbit_list, pole_phase=None):
    orbit_change = 45222 # OrbitChangeNumber
    new_jd0 = 3949.2954391720704734 # orbit change MJD2000

    n_entries = len(orbit_list)

    # combine orbit number and south pole phase, to real number 
    real_orbits = numpy.array(orbit_list)
    if pole_phase != None and len(pole_phase) == n_entries and n_entries > 0:
        real_orbits += pole_phase

    old_jd0           = 790.122714
    old_orbits_per_jd = 501./35.
    new_orbits_per_jd = 431./30.

    old_jds = real_orbits / old_orbits_per_jd + old_jd0
    new_jds = (real_orbits-orbit_change) / new_orbits_per_jd + new_jd0

    offset = date(2000,1,1)
    offset = (offset - date(1,1,1)).days
    jds = numpy.where(real_orbits < orbit_change, old_jds, new_jds) + offset + 3 # +3 what the?
    return jds
