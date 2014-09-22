# -*- coding: iso-8859-1 -*-
#
# COPYRIGHT (c) 2014 SRON (pieter.van.der.meer@sron.nl)
#
#   This is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License, version 2, as
#   published by the Free Software Foundation.
#
#   The software is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, 
#   Boston, MA  02111-1307, USA.
#

"""
Contains functionality related to ENVISAT satellite operations.
"""

import numpy
import time
from datetime import date
import ctypes as ct

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

class PhaseConverter:
    """
    Class that can retrieve polar and eclipse phase information for given julian dates.
    Uses an external library to do its work.
    """

    def __init__(self):
        self.phaselib = ct.CDLL('testlib.so')

    def get_phase(jds, eclipseMode=True):
        """
        input: jds: numpy array of julian dates
               eclipseMode: boolean flag to indicate eclipse mode or polar mode
        returns: numpy array of array phases (same size as 'jds')
        """
        n_records = jds.size
        saa_flag = False;
        dummy_orbit = 0
        # prepare orbit phase output array for c
        phases = numpy.zeros(n_records) # TODO: directly specify dtype=float64 or so?
        phases = phases.astype(numpy.float32)
        phases_c = phases.ctypes.data_as(ct.POINTER(ct.c_float))
        # prepare julianday input array for c
        jds = jds.astype(numpy.float64)
        jds_c = jds.ctypes.data_as(ct.POINTER(ct.c_double))

        phaselib._GET_SCIA_ROE_ORBITPHASE(ct.c_bool(eclipseMode), ct.c_bool(saa_flag), ct.c_long(n_records), ct.c_long(dummy_orbit), phases_c, jds_c) 
        return phases
