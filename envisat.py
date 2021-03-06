#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains functionality related to ENVISAT satellite operations.
"""

from __future__ import print_function, division

import numpy as np
import time
from datetime import date
import ctypes as ct
import re

#-- variables ------------------------------------------------------------------

last_orbit = 53200

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
        if v1 < 1 or v2 > last_orbit:
            raise ArgumentTypeError( msg2 )
        return (v1, v2)
    else:
        return v1

# returns MJD dates (starting 1-1-1) for given ENVISAT absolute orbit numbers
def convert_orbit_to_jd(orbit_list, pole_phase=None):
    orbit_change = 45222 # OrbitChangeNumber
    new_jd0 = 3949.2954391720704734 # orbit change MJD2000

    n_entries = len(orbit_list)

    # combine orbit number and south pole phase, to real number 
    real_orbits = np.array(orbit_list)
    if pole_phase != None and len(pole_phase) == n_entries and n_entries > 0:
        real_orbits += pole_phase

    old_jd0           = 790.122714
    old_orbits_per_jd = 501./35.
    new_orbits_per_jd = 431./30.

    old_jds = real_orbits / old_orbits_per_jd + old_jd0
    new_jds = (real_orbits-orbit_change) / new_orbits_per_jd + new_jd0

    offset = date(2000,1,1)
    offset = (offset - date(1,1,1)).days
    jds = np.where(real_orbits < orbit_change, old_jds, new_jds) + offset + 3 # +3 what the?
    return jds

class PhaseConverter:
    """
    Class that can retrieve polar and eclipse phase information for given julian dates.
    Uses an external library to do its work.
    """

    def __init__(self):
        self.phaselib = ct.CDLL('phaselib.so')

    def get_phase(self, jds, eclipseMode=True, getOrbits=False):
        """
        Convert MJD2000 to accurate orbit number and orbit phase (eclipse-based or pole-based).

        Parameters
        ----------

        jds : numpy array, float64
            julian dates or a float64 scalar
        eclipseMode : bool, optional
            flag to indicate eclipse mode or polar mode
        getOrbits : bool, optional
            also return orbit numbers

        Returns
        -------

        numpy array of array phases (same size as 'jds')

        or

        the phases and an array or absolute orbit numbers
         
        """
        if not isinstance(jds, np.ndarray):
            jds = np.array(jds)

        n_records = jds.size
        saa_flag = False;
        dummy_orbit = 0
        # prepare orbit phase output array for c
        phases = np.zeros(n_records) # TODO: directly specify dtype=float64 or so?
        phases = phases.astype(np.float32)
        phases_c = phases.ctypes.data_as(ct.POINTER(ct.c_float))
        # prepare julianday input array for c
        jds = jds.astype(np.float64)
        jds_c = jds.ctypes.data_as(ct.POINTER(ct.c_double))

        # orbits = np.empty(n_records)
        # orbits_c = orbits.ctypes.data_as(ct.POINTER(ct.c_long))
        orbits = np.empty(n_records) # TODO: directly specify dtype=float64 or so?
        orbits = orbits.astype(np.int32)
        orbits_c = orbits.ctypes.data_as(ct.POINTER(ct.c_long))

        if getOrbits:
            self.phaselib._GET_SCIA_ROE_ORBITPHASE_ORBIT(ct.c_bool(eclipseMode), ct.c_bool(saa_flag), ct.c_long(n_records), orbits_c, phases_c, jds_c) 
            return np.float64(phases), orbits
        else:
            self.phaselib._GET_SCIA_ROE_ORBITPHASE(ct.c_bool(eclipseMode), ct.c_bool(saa_flag), ct.c_long(n_records), ct.c_long(dummy_orbit), phases_c, jds_c) 
            return phases

