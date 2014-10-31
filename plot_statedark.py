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

from __future__ import print_function
from __future__ import division

from read_statedark_module import sdmf_read_statedark
import numpy
import h5py
from numpy import arange, sin, pi
import matplotlib
from matplotlib.figure import Figure
from pylab import plot, savefig

db_name = "/Users/pieterm/sdmf_statedark_ch8.h5"
orbit_range = [24000,24000]
pix = 140 

n_orbits = orbit_range[1] - orbit_range[0]

data = sdmf_read_statedark(orbit_range, 26, db_name, darksignal=True, darkmedian=True, stddev=True)

#print(data)

#print(data['DarkSignal'])
#print(data['DarkMedian'])
#print(data['StdDev'])

# orbits, pixels, states
#print((data['DarkSignal'])[:,140])
print(data['DarkSignal'].shape)

#hor_axis = numpy.arange(orbit_range[0], orbit_range[1])
#idx = numpy.nonzero(vals)
#print(vals)

colors = ['b','g','r','c','m','y','k','w']

for i in range((data['DarkSignal']).shape[2]):
    vals = (data['DarkSignal'])[:,pix,i]
    print(vals)
    hor_axis = numpy.arange(vals.shape[0])
    #print(hor_axis.shape)
    plot(hor_axis, vals, color=colors[i%8])

savefig('test.png', dpi=100)
