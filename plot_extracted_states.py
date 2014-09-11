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

import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy
from numpy import arange, sin, pi
import h5py

# Used to guarantee to use at least Wx2.8
import wxversion
#wxversion.ensureMinimal('2.8')
# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')


# reads extracted state executions from database
def read_extracted_states(orbitrange, state_id, calib_db, readoutMean=False, readoutNoise=False, orbitList=False):

    if len(orbitrange) is not 2:
        print('read_extracted_states: orbitrange should have 2 elements')
        return

    dict = {}

    #
    # obtain indices to entries of the requested orbit and state
    #

    fid = h5py.File(calib_db, 'r') # generates an exception
    gid = fid["State_"+str('%02d'%state_id)] # generates an exception
    mt_did = gid['metaTable']
    orbitlist = (gid['orbitList'])[:]
    print('orbitrange=', orbitrange)
#    metaindx = numpy.where((orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1]))
    metaindx = (orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1])
    #print(metaindx)
    if metaindx[0].size is 0:
        print('read_extracted_states: orbit range not present in database')
        dict['status'] = -1
        return dict

    mtbl = mt_did[metaindx]
    dict['mtbl'] = mtbl

    if readoutNoise:
        ds_did      = gid['readoutMean']
        dict['readoutMean'] = ds_did[metaindx,:]

    if readoutNoise:
        ds_did       = gid['readoutNoise']
        dict['readoutNoise'] = ds_did[metaindx,:]

    if orbitList:
        ds_did = gid['orbitList']
        dict['orbitList'] = ds_did[metaindx]

    fid.close()
    dict['status'] = 0
    return(dict)


db_name = "/SCIA/SDMF31/sdmf_extract_calib.h5"
orbit_range = [24000,25000]
pix = 7*1024+140 

n_orbits = orbit_range[1] - orbit_range[0]

data = read_extracted_states(orbit_range, 26, db_name, readoutMean=True, readoutNoise=True, orbitList=True)

#print(data)

# orbits, pixels, states
print(data['readoutMean'].shape)

#hor_axis = numpy.arange(orbit_range[0], orbit_range[1])
#idx = numpy.nonzero(vals)
#print(vals)

colors = ['b','g','r','c','m','y','k','w']

#print(data['mtbl'])
orbits = data['orbitList']
print(orbits)
uorbits = numpy.unique(orbits)
readoutMean = data['readoutMean']
mtbl = data['mtbl']
phases = mtbl['orbitPhase'][:]

for i in range(uorbits.size):
    #print(i, uorbits[i])
    idx = numpy.where(orbits == uorbits[i])
    vals = readoutMean[idx, pix]
    #print(vals)
    vals = vals.flatten()
    hor_axis = phases[idx] #numpy.arange(vals.size)
    #print(vals.shape)
    #print(vals)
    #print(hor_axis.shape)
    #print(hor_axis)
    plt.scatter(hor_axis, vals, color=colors[i%8])

#plt.savefig('test.png', dpi=100)
plt.show()
