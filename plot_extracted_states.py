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

import argparse
import matplotlib
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy
from numpy import arange, sin, pi
import h5py
import ctypes as ct
from sciamachy_module import NonlinCorrector

# Used to guarantee to use at least Wx2.8
import wxversion
#wxversion.ensureMinimal('2.8')
# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')

def onclick(event):
    print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (event.button, event.x, event.y, event.xdata, event.ydata))

# read simudark data for specified orbit (only ch8 data!)
def read_simudark(orbit, ao=None, lc=None, amp1=None):
    dict = {}

    #
    # obtain indices to entries of the requested orbit
    #

    fid = h5py.File("/SCIA/share/SDMF/3.0/sdmf_simudark.h5", 'r') # generates an exception
    gid = fid["ch8"] # generates an exception
    mt_did = gid['metaTable']
    orbitlist = (gid['orbitList'])[:]
    if orbit not in orbitlist:
        print('read_simudark(): orbit not present in database')
        dict['status'] = -1
        return dict

    metaindx = numpy.where(orbitlist == orbit)
    metaindx = numpy.asscalar(numpy.array(metaindx))
    print(metaindx)

    mt_did = gid['metaTable']
    mtbl = mt_did[metaindx]
    dict['mtbl'] = mtbl

    if ao:
        ds_did = gid['ao']
        dict['ao'] = ds_did[:,metaindx]

    if lc:
        ds_did = gid['lc']
        dict['lc'] = ds_did[:,metaindx]

    if amp1:
        ds_did = gid['amp1']
        dict['amp1'] = ds_did[:,metaindx]

    dict['status'] = 0
    return(dict)    

# reads extracted state executions from database
def read_extracted_states(orbitrange, state_id, calib_db, in_orbitlist=None, readoutMean=False, readoutNoise=False, orbitList=False):

    if in_orbitlist is None:
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
    if in_orbitlist is not None:
        #metaindx = orbitlist = in_orbitlist
        metaindx = numpy.in1d(orbitlist, in_orbitlist)
    else:
        print('orbitrange=', orbitrange)
        metaindx = (orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1])
        print(metaindx)

    if metaindx[0].size is 0:
        print('read_extracted_states: orbit range not present in database')
        dict['status'] = -1
        return dict

    mtbl = mt_did[metaindx]
    dict['mtbl'] = mtbl

    if readoutMean:
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

#-------------------------------------------------------------------------------
# computes the variable part of the simudark product for specified orbital 
# phases.
#
# INPUT: (structure)
#       - phases:    orbital phases [0.0..1.0] (1D array)
#       - pet:       pixel exposure time
#       - amp1:      amplitude of fundamental frequency (1D array)
#       - amp2:      relative amplitude of 1st harmonic
#       - phase1:    phase shift of fundamental frequency
#       - phase2:    phase shift of 1st harmonic
#
# KEYWORDS:
#       - n_harmonics: input: how many harmonics to use (0 or 1, default 1)
#
# RETURNS:
#       - variable part per pixel and phase (2D array)
#-------------------------------------------------------------------------------
def simudark_orbvar_function(d, n_harmonics=1):

    exec_count = d['phases'].size
    n_pixels   = d['amp1'].size
    func       = numpy.empty((n_pixels, exec_count))

    pet_amp   = numpy.matrix(d['pet'] * d['amp1'])
    amp2      = d['amp2'] * pet_amp
    # TODO: expand amp2 to channel width if not already so?!
    #if n_elements(amp2) eq n_pixels/1024 then $
    #    amp2 = rebin(amp2,n_pixels,/sample)

    # correct fundamental frequency
    phases1 = numpy.matrix(d['phases']+d['phase1']) # orbital phases shifted by fundamental phase shift
    func = pet_amp.T * numpy.cos(2*numpy.pi*phases1)

    # correct 1st harmonic
    if n_harmonics >= 1:
        phases2 = numpy.matrix(d['phases']+d['phase2'])
        func += amp2.T * numpy.cos(4*numpy.pi*phases2)

    return func.T

# just a smoke test.. 
def test_simudark_orbvar_function():
    d = {}
    d['phases'] = numpy.array([0.1,0.2,0.3]) # 3 state executions
    d['pet'] = 0.998
    d['amp1'] = numpy.array([1,1.1,1,1.1,1,1.1,1,1.1,1,1.1]) # 10 pixels
    d['amp2'] = .12
    d['phase1'] = .1
    d['phase2'] = .2
    print(simudark_orbvar_function(d))

# test function to check old simudark quality
# pixnr: pixel in channel 8 [0..1023]
def check_eclipse_calib(pixnr):
    orbit = 24000

    #
    # get state executions and correct them for nonlinearity
    # 

    db_name = "/SCIA/SDMF31/sdmf_extract_calib.h5"
    states = read_extracted_states([orbit,orbit], 8, db_name, readoutMean=True)
    state_mtbl = states['mtbl']
    print(states['readoutMean'].shape)
    readouts = states['readoutMean'] #states['readoutMean'][:,7*1024:8*1024]
    print(readouts.shape)
    nlc = NonlinCorrector()
    for idx_exec in range(readouts.shape[0]):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])

    #
    # correct for dark current with simudark v1
    #

    petcorr = 1.18125e-3
    simudark = read_simudark(orbit, ao=True, lc=True, amp1=True)
    simudark_mtbl = simudark['mtbl']
    state_phases = state_mtbl['orbitPhase'][:]
    readouts_ch8 = readouts[:,7*1024:8*1024]
    readouts_ch8 -= simudark['ao']
    readouts_ch8 -= simudark['lc'] * (1-petcorr) # TODO: PET belonging to the state data
    d = {}
    d['phases'] = state_phases
    d['pet'] = numpy.array(1.0-petcorr) * numpy.ones(1024) # TODO: get pet from state data!
    d['amp1'] = simudark['amp1'][:]
    d['amp2'] = simudark_mtbl['AMP2'] #[:]
    d['phase1'] = simudark_mtbl['PHASE1'] #[:]
    d['phase2'] = simudark_mtbl['PHASE2'] #[:]
    funk = numpy.array(simudark_orbvar_function(d))
    print(funk.shape)
    print(d)
    #plt.scatter(state_phases, readouts_ch8[:,pixnr]-funk[:,pixnr])
    fig = plt.figure()
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    darklevel = (simudark['ao']+simudark['lc']*(1.0-petcorr))[pixnr]
    plt.suptitle("dark level = "+str(darklevel))
    plt.scatter(state_phases, readouts_ch8[:,pixnr], color='g')
    plt.scatter(state_phases, funk[:,pixnr], color='r')
    plt.show()
    return

# plot dark state executions (detrend and overplot them)
def check_darkstates():
    testlib = ct.CDLL('testlib.so')

    #testlib.myprint()
    #d = ct.c_double(1.234)
    #ds = (ct.c_double * 4)()
    #for i in range(4):
    #    ds[i] = i
    #testlib.print_double(d)
    #testlib.print_doubles(ds, ct.c_long(len(ds)))

    #n_elem = 4
    #zs = numpy.arange(n_elem)
    #print(zs)
    #zs = zs.astype(numpy.float64)
    #testlib.print_doubles(zs.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_long(n_elem))

    db_name = "/SCIA/SDMF31/sdmf_extract_calib.h5"
    orbit_range = [24000,24005]
    pix = 7*1024+140 

    n_orbits = orbit_range[1] - orbit_range[0]

    in_orbitlist =  [20237,20238,20666,20667,21081,21082,21511,21512,21954,21955,22384,22385,22799,22800,23214,23215,23629,23630,24044,24045,24474,24475,24889,24890,25318,25319,25733,25734,26163,26164,26563,26564,27022,27023,27437,27438,28052,28053,28295,28296]

    #data = read_extracted_states(orbit_range, 8, db_name, readoutMean=True, readoutNoise=True, orbitList=True)
    data = read_extracted_states(orbit_range, 8, db_name, in_orbitlist=in_orbitlist, readoutMean=True, readoutNoise=True, orbitList=True)

    #print(data)

    # orbits, pixels, states
    print(data['readoutMean'].shape)

    #hor_axis = numpy.arange(orbit_range[0], orbit_range[1])
    #idx = numpy.nonzero(vals)
    #print(vals)

    colors = ['b','g','r','c','m','y','k']

    #print(data['mtbl'])
    orbits = data['orbitList']
    print(orbits)
    uorbits = numpy.unique(orbits)
    readoutMean = data['readoutMean']
    mtbl = data['mtbl']
    n_records = mtbl.size
    jds = mtbl['julianDay'][:]

    #
    # get eclipse orbit phases instead of the SDMF 3.1 polar-based ones
    #

    saa_flag = False;
    dummy_orbit = 0
    # prepare orbit phase output array for c
    e_phases = numpy.zeros(n_records) # TODO: directly specify dtype=float64 or so?
    e_phases = e_phases.astype(numpy.float32)
    e_phases_c = e_phases.ctypes.data_as(ct.POINTER(ct.c_float))
    # prepare julianday input array for c
    jds = jds.astype(numpy.float64)
    jds_c = jds.ctypes.data_as(ct.POINTER(ct.c_double))
    testlib._GET_SCIA_ROE_ORBITPHASE(ct.c_bool(True), ct.c_bool(saa_flag), ct.c_long(n_records), ct.c_long(dummy_orbit), e_phases_c, jds_c) 

    print(e_phases)

    phases = mtbl['orbitPhase'][:]

    for i in range(uorbits.size):
        #print(i, uorbits[i])
        idx = numpy.where(orbits == uorbits[i])
        vals = readoutMean[idx, pix]
        #print(vals)
        vals = vals.flatten()
        #hor_axis = phases[idx] #numpy.arange(vals.size)
        hor_axis = e_phases[idx] #numpy.arange(vals.size)
        #print(vals.shape)
        #print(vals)
        #print(hor_axis.shape)
        #print(hor_axis)
        plt.scatter(hor_axis, vals-numpy.mean(vals), color=colors[i%7])
        #plt.scatter(hor_axis, vals, color=colors[i%7])
     
    #plt.savefig('test.png', dpi=100)
    plt.show()
    return

#
# main
#

parser = argparse.ArgumentParser()
parser.add_argument("pixnr")
args = parser.parse_args()
print(args.pixnr)

#test_simudark_orbvar_function()
check_eclipse_calib(int(args.pixnr))
#check_darkstates()

