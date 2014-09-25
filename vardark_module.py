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

import h5py
import numpy
import numpy as np
from numpy import cos, pi
from kapteyn import kmpfit
import matplotlib.pyplot as plt

from envisat import PhaseConverter
from sciamachy_module import NonlinCorrector, read_extracted_states, petcorr, orbitfilter

#-------------------------SECTION VERSION-----------------------------------

_swVersion = {'major': 0,
              'minor': 2,
              'revision' : 0}
_calibVersion = {'major': 0,
                 'minor': 1,
                 'revision' : 0}
_dbVersion = {'major': 0,
              'minor': 1,
              'revision' : 0}

#- globals ---------------------------------------------------------------------

nlc = NonlinCorrector()
phaseconv = PhaseConverter()
ofilt = orbitfilter()
fname = '/SCIA/SDMF31/sdmf_extract_calib.h5'
n_pix = 1024

#- functions -------------------------------------------------------------------

def scia_dark_fun1(p, x):
    ao = p[0]
    dc = p[1]
    amp1 = p[2] # wave amp, relative to dark
    trend = p[3] # relative to dark
    phase_shift1 = p[4]
    #amp2 = p[5]
    #phase_shift2 = p[6]
    #print(ao,dc,amp1,trend,phase_shift1)
    # Extract exposure information and orbit phase from x
    orbit_phase, pet, coadd = x
    #print(orbit_phase)
    n_x = orbit_phase.size # nr of datapoints

    dark = numpy.zeros(n_x)
    dark += dc
    #print(amp1 * cos(2*pi*(orbit_phase + phase_shift1)))
    dark += amp1 * cos(2*pi*(orbit_phase + phase_shift1))
    #print(amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2)))
    #dark += amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2))
    dark += trend * orbit_phase

#    return (dark*pet + ao)*coadd    
    return dark*pet + ao

def scia_dark_residuals1(p, data):
    x, y, yerr = data 
    return y - scia_dark_fun1(p, x)

def scia_dark_residuals1e(p, data):
    x, y, yerr = data 
    return (y - scia_dark_fun1(p, x)) / yerr

def extract_two_dark_states(orbit, stateid1, stateid2):
    # orbit 1
    orbrange = [orbit,orbit]
    states = read_extracted_states(orbrange, stateid1, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states['mtbl']
    jds1_8 = state_mtbl['julianDay'][:]
    #print(states['readoutMean'].shape)
    readouts = states['readoutMean']
    noise = states['readoutNoise']
    n_exec1_8 = readouts.shape[0]
    for idx_exec in range(n_exec1_8):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    pet8 = states['pet'][0]
    #print('pet8=',pet8)
    coadd8 = states['coadd'][0]
    readouts8_ = readouts[:,7*1024:8*1024] #/ coadd (was already done)
    noise8_ = noise[:,7*1024:8*1024] / numpy.sqrt(coadd8)
    state8_phases = state_mtbl['orbitPhase'][:]

    states63 = read_extracted_states(orbrange, stateid2, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states63['mtbl']
    jds1_63 = state_mtbl['julianDay'][:]
    #print(states['readoutMean'].shape)
    readouts = states63['readoutMean']
    noise = states63['readoutNoise']
    n_exec1_63 = readouts.shape[0]
    for idx_exec in range(n_exec1_63):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    pet63 = states63['pet'][0]
    #print('pet63=',pet63)
    coadd63 = states['coadd'][0]
    readouts63_ = readouts[:,7*1024:8*1024] #/ coadd (was already done)
    noise63_ = noise[:,7*1024:8*1024] / numpy.sqrt(coadd63)
    state63_phases = state_mtbl['orbitPhase'][:]

    # glue data from first orbit
    n_exec1 = n_exec1_8 + n_exec1_63
    state_phases1 = numpy.concatenate((state8_phases, state63_phases)) 
    readouts1 = numpy.concatenate((readouts8_, readouts63_))
    noise1 = numpy.concatenate((noise8_, noise63_))

    # orbit 2
    orbrange = [orbit+1,orbit+1]
    states = read_extracted_states(orbrange, stateid1, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states['mtbl']
    jds2_8 = state_mtbl['julianDay'][:]
    #print(states['readoutMean'].shape)
    readouts = states['readoutMean']
    noise = states['readoutNoise']
    n_exec2_8 = readouts.shape[0]
    for idx_exec in range(n_exec2_8):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    readouts2_8 = readouts[:,7*1024:8*1024] #/ coadd (was already done)
    noise2_8 = noise[:,7*1024:8*1024] / numpy.sqrt(coadd8)
    state_phases2 = state_mtbl['orbitPhase'][:]

    states63 = read_extracted_states(orbrange, stateid2, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states63['mtbl']
    jds2_63 = state_mtbl['julianDay'][:]
    #print(states['readoutMean'].shape)
    readouts = states63['readoutMean']
    noise = states63['readoutNoise']
    n_exec2_63 = readouts.shape[0]
    for idx_exec in range(n_exec2_63):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    pet63 = states63['pet'][0]
    #print('pet63=',pet63)
    coadd63 = states['coadd'][0]
    readouts63_ = readouts[:,7*1024:8*1024] #/ coadd (was already done)
    noise63_ = noise[:,7*1024:8*1024] / numpy.sqrt(coadd63)
    state63_phases = state_mtbl['orbitPhase'][:]

    # glue data from second orbit
    n_exec2 = n_exec2_8 + n_exec2_63
    state_phases2 = numpy.concatenate((state_phases2, state63_phases)) 
    readouts2 = numpy.concatenate((readouts2_8, readouts63_))
    noise2 = numpy.concatenate((noise2_8, noise63_))

    # combine orbit data
    n_exec = n_exec1 + n_exec2
    all_state_phases = numpy.concatenate((state_phases1, state_phases2+1.0))
    all_readouts = numpy.concatenate((readouts1, readouts2))
    all_sigmas = numpy.concatenate((noise1, noise2))

    pet = numpy.concatenate((numpy.zeros(n_exec1_8)+pet8, numpy.zeros(n_exec1_63)+pet63, numpy.zeros(n_exec2_8)+pet8, numpy.zeros(n_exec2_63)+pet63))
    coadd = numpy.concatenate((numpy.zeros(n_exec1_8)+coadd8, numpy.zeros(n_exec1_63)+coadd63, numpy.zeros(n_exec2_8)+coadd8, numpy.zeros(n_exec2_63)+coadd63))
    jds = numpy.concatenate((jds1_8, jds1_63, jds2_8, jds2_63))

    # convert all juliandates to eclipse phases
    #print(jds.size, jds)
    ephases = phaseconv.get_phase(jds)
    ephases[n_exec1:n_exec] += 1.0 # second orbit should have orbit phase +1

    return n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases

# TODO: later orbits have different dark state definitions..
def extract_05_10_dark_states(orbit):
    return extract_two_dark_states(orbit, 8, 63)

# TODO: later orbits have different dark state definitions..
def extract_short_dark_states(orbit):
    return extract_two_dark_states(orbit, 26, 46)

def extract_dark_states(orbit, shortFlag=False, longFlag=True):
    #print(shortFlag, longFlag)
    if shortFlag:
        n_exec_s, all_state_phases_s, pet_s, coadd_s, all_readouts_s, all_sigmas_s, ephases_s = extract_short_dark_states(orbit)
    if longFlag:
        n_exec_l, all_state_phases_l, pet_l, coadd_l, all_readouts_l, all_sigmas_l, ephases_l = extract_05_10_dark_states(orbit)

    if shortFlag and longFlag:
        n_exec = n_exec_s+n_exec_l
        all_state_phases = numpy.concatenate((all_state_phases_s, all_state_phases_l))
        pet = numpy.concatenate((pet_s, pet_l))
        coadd = numpy.concatenate((coadd_s, coadd_l))
        all_readouts = numpy.concatenate((all_readouts_s, all_readouts_l))
        all_sigmas = numpy.concatenate((all_sigmas_s, all_sigmas_l))
        ephases = numpy.concatenate((ephases_s, ephases_l))
    elif shortFlag and not longFlag:
        n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = n_exec_s, all_state_phases_s, pet_s, coadd_s, all_readouts_s, all_sigmas_s, ephases_s
    elif not shortFlag and longFlag:
        n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = n_exec_l, all_state_phases_l, pet_l, coadd_l, all_readouts_l, all_sigmas_l, ephases_l

    return n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases

# fit dark model to two neighbouring monthly calibration orbits
#def fit_monthly(orbit, shortFlag=False, longFlag=True):
def fit_monthly(orbit, verbose=False, **kwargs):

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = extract_dark_states(orbit, **kwargs)

    #
    # fit it
    #

    x = all_state_phases, pet, coadd
    if verbose:
        print(all_state_phases.shape)
        print(pet.shape)
        print(pet)
        print(coadd)
    #fitobj = kmpfit.Fitter(residuals=scia_dark_residuals1, data=(x, all_readouts[:,pixnr], all_sigmas[:,pixnr]))
    
    aoinfo = dict(fixed=False, limits=[0,10000])
    lcinfo = dict(fixed=False, limits=[-100000.,+100000.])
    amp1info = dict(fixed=False, limits=[-1000,+1000])
    trendinfo = dict(fixed=False, limits=[-1000,+1000])
    amp2info = dict(fixed=False, limits=[-1.,+1.])
    phase1info = dict(fixed=False, limits=[-3.,+3.])
    phase2info = dict(fixed=False, limits=[-3.,+3.])
    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info] 

    # prepare initial parameters
    ao0 = 4000
    dc0 = 4000
    amp1_0 = 10
    trend0 = 0
    amp2_0 = 0.1
    phase_offset1 = 0
    phase_offset2 = 0
    p0 = numpy.array([ao0, dc0, amp1_0, trend0, phase_offset1]) #, phase_offset2

    # ..and fit
    n_done = 0
    statuses = numpy.zeros(n_pix)
    res_phases = numpy.zeros(n_pix)
    for pixnr in range(n_pix):
        #print(pixnr)
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        #print('readouts = ', pix_readouts)
        #print('sigmas = ', pix_sigmas)
        #print(numpy.sum(pix_sigmas))
        if not numpy.isnan(numpy.sum(pix_readouts)) and numpy.all(pix_sigmas != 0):
            fitobj = kmpfit.simplefit(scia_dark_fun1, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            n_done += 1
        else:
            continue
        #fitobj.fit(params0=p0)
        if verbose:
            print(pixnr, fitobj.params[4] % 1, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_phases[pixnr] = fitobj.params[4]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()
        #else:
        #   print("Optimal parameters: ", fitobj.params)

    channel_phase = numpy.median(res_phases[numpy.where(statuses >= 0) and numpy.where(res_phases != -3)] % 1)
    if verbose:
        print('channel median phase =', channel_phase)
    phase1info = dict(fixed=True, limits=[-3.,+3.])
    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info] 
    p0[4] = channel_phase

    #
    # pass 2 - fix phase shift
    #

    aos = numpy.zeros(n_pix)
    lcs = numpy.zeros(n_pix)
    amps = numpy.zeros(n_pix)
    trends = numpy.zeros(n_pix)
    for pixnr in range(n_pix):
        #print(pixnr)
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        #print('readouts = ', pix_readouts)
        #print('sigmas = ', pix_sigmas)
        #print(numpy.sum(pix_sigmas))
        if not numpy.isnan(numpy.sum(pix_readouts)) and numpy.all(pix_sigmas != 0):
            fitobj = kmpfit.simplefit(scia_dark_fun1, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            n_done += 1
        else:
            continue
        #fitobj.fit(params0=p0)
        if verbose:
            print(pixnr, fitobj.message)
        statuses[pixnr] = fitobj.status
        aos[pixnr] = fitobj.params[0]
        lcs[pixnr] = fitobj.params[1]
        amps[pixnr] = fitobj.params[2]
        trends[pixnr] = fitobj.params[3]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()

    #plt.cla()
    #plt.scatter(all_state_phases, all_readouts[:,pixnr])
    #plt.show()

    return channel_phase, aos, lcs, amps, trends

# fit dark model to an orbits with normal (eclipse) dark states
# use parameters computed from nearest calibration orbits
def fit_eclipse_orbit(orbit, aos, lcs, amps, channel_phaseshift, verbose=False, **kwargs):

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = extract_dark_states(orbit, **kwargs)

    #
    # fit it
    #

    x = all_state_phases, pet, coadd
    #print(all_state_phases.shape)
    #print(pet.shape)
    #print(pet)
    #print(coadd.shape)
    #fitobj = kmpfit.Fitter(residuals=scia_dark_residuals1, data=(x, all_readouts[:,pixnr], all_sigmas[:,pixnr]))
    
    aoinfo = dict(fixed=True, limits=[0,10000])
    lcinfo = dict(fixed=False, limits=[-100000.,+100000.])
    amp1info = dict(fixed=True, limits=[-1000,+1000])
    trendinfo = dict(fixed=False, limits=[-1000,+1000])
    amp2info = dict(fixed=False, limits=[-1.,+1.])
    phase1info = dict(fixed=True, limits=[-3.,+3.])
    phase2info = dict(fixed=False, limits=[-3.,+3.])
    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info] 

    # ..and fit
    n_done = 0
    statuses = numpy.zeros(n_pix)
    res_phases = numpy.zeros(n_pix) + numpy.nan
    res_trends = numpy.zeros(n_pix) + numpy.nan
    for pixnr in range(n_pix):
        #print(pixnr)
        # prepare initial parameters
        p0 = numpy.array([aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phaseshift]) 

        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        #print('readouts = ', pix_readouts)
        #print('sigmas = ', pix_sigmas)
        #print(numpy.sum(pix_sigmas))
        if (aos[pixnr] is not 0) and (not numpy.isnan(numpy.sum(pix_readouts))) and (numpy.all(pix_sigmas != 0)):
            fitobj = kmpfit.simplefit(scia_dark_fun1, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            n_done += 1
        else:
            continue
        #fitobj.fit(params0=p0)
        if verbose:
            print(pixnr, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_trends[pixnr] = fitobj.params[3]
        res_phases[pixnr] = fitobj.params[4]
        if (fitobj.status <= 0):
           #print('Error message = ', fitobj.message)
           quit()
        #else:
        #   print("Optimal parameters: ", fitobj.params)

    print("done")

    return x, lcs, res_trends, all_readouts, all_sigmas

# compute thermal background trend between two normal orbits 
# also computes actual thermal background offset (excluding trend or oscillation)
def compute_trend(orbit, aos, amps, phaseshift, **kwargs):

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = extract_dark_states(orbit, **kwargs)

    # divide by coadding factor, subtract analog offset and divide by exposure time to get thermal background in BU/s
    thermal_background = all_readouts
    #thermal_background /= numpy.matrix(coadd).T * numpy.ones(n_pix) # this was already done.. 
    thermal_background -= (numpy.matrix(aos).T * numpy.ones(n_exec)).T
    thermal_background /= numpy.matrix(pet).T * numpy.ones(n_pix)

    trends = numpy.empty(n_pix) + numpy.nan
    lcs = numpy.empty(n_pix) + numpy.nan

    # loop over all channel pixels
    for pixnr in range(n_pix):
        if aos[pixnr] is 0:
            continue
        if numpy.all(all_sigmas[:,pixnr] is 0):
            continue 
        background = thermal_background[:,pixnr]
        if numpy.isnan(numpy.sum(background)):
            continue 
        idx1 = numpy.where(ephases < 1.0)
        idx2 = numpy.where(ephases >= 1.0)
        avg1 = numpy.mean(background[idx1])
        avg2 = numpy.mean(background[idx2])
        phi1 = numpy.mean(ephases[idx1])
        phi2 = numpy.mean(ephases[idx2])
        lcs[pixnr] = avg1 - amps[pixnr]*cos(2*pi*(phi1+phaseshift))

        trends[pixnr] = (avg2-avg1)/(phi2-phi1)

    return trends, lcs

class VarDark:
    def __init__(self, orbit,ao, lc, amp, trend, phaseshift):
        self.numPixels = 1024
        self.absOrbit = orbit
        self.ao = ao
        self.lc = lc
        self.amp = amp
        self.trend = trend
        self.phaseshift = phaseshift
        self.mtbl = dict()
        self.mtbl['julianDay'] = 1.234 # TODO convert orbit to JD
        self.mtbl['obmTemp'] = 999 # TODO
        self.mtbl['detTemp'] = numpy.array([999,999,999,999,999,999,999,999]) # TODO

class VarDarkDb:
    def __init__( self, args=None, db_name='./sdmf_vardark.h5',
                  truncate=False, calibration=None, verbose=False ):
        self.monthly_orbit = -1
        self.created = False
        if args:
            self.db_name  = args.db_name
            self.calibration = args.calibration.copy()
            self.truncate = args.truncate
            self.verbose  = args.verbose
        else:
            self.db_name  = db_name
            self.calibration = calibration
            self.truncate = truncate
            self.verbose  = verbose

    def checkDataBase(self):
        with h5py.File( self.db_name, 'r' ) as fid:
            mystr = ','.join(list(self.calibration.astype('str')))
            if fid.attrs['calibOptions'] != mystr:
                print( 'Fatal:', 'incompatible calibration options' )
                raise dbError('incompatible calibration options')
            myversion = '%(major)d.%(minor)d' % _swVersion
            if fid.attrs['swVersion'].rsplit('.',1)[0] != myversion:
                print( 'Fatal:', 'incompatible with _swVersion' )
                raise dbError('incompatible with _swVersion')
            myversion = '%(major)d.%(minor)d' % _dbVersion
            if fid.attrs['dbVersion'].rsplit('.',1)[0] != myversion:
                print( 'Fatal:', 'incompatible with _dbVersion' )
                raise dbError('incompatible with _dbVersion')
            myversion = '%(major)d.%(minor)d' % _calibVersion
            if fid.attrs['calibVersion'].rsplit('.',1)[0] != myversion:
                print( 'Fatal:', 'incompatible with _calibVersion' )
                raise dbError('incompatible with _calibVersion')

    def fill_mtbl(self, vardark):
        from datetime import datetime 

        fmtMTBL  = 'float64,a20,uint16,uint16,float32,float32,8float32'
        nameMTBL = ('julianDay','entryDate','absOrbit','quality','phaseShift','obmTemp','detTemp')

        self.mtbl = np.empty(1, dtype=fmtMTBL )
        self.mtbl.dtype.names = nameMTBL
        self.mtbl['julianDay'] = vardark.mtbl['julianDay']
        self.mtbl['entryDate'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.mtbl['absOrbit'] = vardark.absOrbit
        self.mtbl['quality'] = 0
        self.mtbl['phaseShift'] = vardark.phaseshift
        self.mtbl['obmTemp'] = vardark.mtbl['obmTemp']
        self.mtbl['detTemp'] = vardark.mtbl['detTemp']

    def create(self, vardark):
        with h5py.File( self.db_name, 'w', libver='latest' ) as fid:
            n_pix = vardark.numPixels
            if isinstance(vardark.absOrbit, numpy.ndarray):
                reshaped_orbits = vardark.absOrbit.reshape(1,)
            else:
                reshaped_orbits = numpy.array([vardark.absOrbit])
            ds = fid.create_dataset( 'orbitList', dtype='uint16',
                                     data=reshaped_orbits, 
                                     maxshape=(None,), chunks=(512,) )
            ds = fid.create_dataset( 'metaTable', 
                                     data=self.mtbl,
                                     chunks=(16384 // self.mtbl.dtype.itemsize,),
                                     shuffle=True, compression='gzip',
                                     compression_opts=1, maxshape=(None,) )
            ds = fid.create_dataset( 'ao', 
                                     data=vardark.ao.reshape(1,n_pix), 
                                     maxshape=(None,n_pix), 
                                     chunks=(8,n_pix), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'lc', 
                                     data=vardark.lc.reshape(1,n_pix), 
                                     maxshape=(None,n_pix), 
                                     chunks=(8,n_pix), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'amp', 
                                     data=vardark.amp.reshape(1,n_pix), 
                                     maxshape=(None,n_pix), 
                                     chunks=(8,n_pix), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'trend', 
                                     data=vardark.trend.reshape(1,n_pix), 
                                     maxshape=(None,n_pix), 
                                     chunks=(8,n_pix), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )

            # create attributes in the HDF5 root
            if self.calibration is not None:
                mystr = ','.join(list(self.calibration.astype('str')))
                fid.attrs['calibOptions'] = mystr
            fid.attrs['swVersion'] = \
                '%(major)d.%(minor)d.%(revision)d' % _swVersion
            fid.attrs['dbVersion'] = \
                '%(major)d.%(minor)d.%(revision)d' % _dbVersion
            fid.attrs['calibVersion'] = \
                '%(major)d.%(minor)d.%(revision)d' % _calibVersion

            self.created = True

    def append(self, vardark):
        with h5py.File( self.db_name, 'r+' ) as fid:
            dset = fid['orbitList']         # orbitList
            ax1 = dset.len()
            dset.resize(ax1+1, axis=0)
            dset[ax1] = vardark.absOrbit
            orbitList = dset[:]
            dset = fid['metaTable']         # metaTable
            dset.resize(ax1+1, axis=0)
            dset[ax1] = self.mtbl
            dset = fid['ao']               # analog offset
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = vardark.ao.reshape(1,n_pix)
            dset = fid['lc']       # leakage current (constant part of thermal background in BU/s)
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = vardark.lc.reshape(1,n_pix)
            dset = fid['amp']       # orbital variation amplitude
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = vardark.amp.reshape(1,n_pix)
            dset = fid['trend']          # orbit-to-orbit trend (slope)
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = vardark.trend.reshape(1,n_pix)

    def rewrite(self, vardark):
        with h5py.File( self.db_name, 'r+' ) as fid:
            dset = fid['orbitList']         # orbitList
            orbitList = dset[:]
            ax1 = np.nonzero(orbitList == vardark.absOrbit)[0][0]
            dset = fid['metaTable']         # metaTable
            dset[ax1] = self.mtbl
            dset = fid['ao']               # analog offset
            dset[ax1,:] = vardark.lc.reshape(1,n_pix)
            dset = fid['lc']       # leakage current (constant part of thermal background in BU/s)
            dset[ax1,:] = vardark.ao.reshape(1,n_pix)
            dset = fid['amp']       # orbital variation amplitude
            dset[ax1,:] = vardark.amp.reshape(1,n_pix)
            dset = fid['trend']          # orbit-to-orbit trend (slope)
            dset[ax1,:] = vardark.trend.reshape(1,n_pix)

    def store(self, vardark, verbose=False):
        self.fill_mtbl( vardark )
        if not h5py.is_hdf5( self.db_name ):
            if verbose: print( '* Info: create new database' )
            self.create( vardark )
        elif self.truncate: 
            if verbose: print( '* Info: replace database (and start a new)' )
            self.create( vardark )
        else:
            self.checkDataBase()
            with h5py.File( self.db_name, 'r' ) as fid:
                dset = fid['orbitList']
                orbitList = dset[:]

            if np.nonzero(orbitList == vardark.absOrbit)[0].size == 0:
                if verbose: print( '* Info: append new SMR to database' )
                self.append( vardark )
            else:
                if verbose: print( '* Info: overwrite entry in database' )
                self.rewrite( vardark )

    def calc_and_store_orbit(self, orbit, verbose=False, **kwargs):
        normal_orbit = orbit
        monthly_orbit = ofilt.get_closest_monthly(orbit)
        if self.monthly_orbit != monthly_orbit:
            # calculate if monthly if not buffered
            channel_phase, aos, lcs, amps, trends = fit_monthly(monthly_orbit, **kwargs)
            self.monthly_orbit = monthly_orbit
            self.aos = aos
            self.lcs = lcs
            self.amps = amps
            self.trends = trends
            self.channel_phase = channel_phase 
        else:
            # just use the buffered version
            aos = self.aos
            lcs = self.lcs
            amps = self.amps
            trends = self.trends
            channel_phase = self.channel_phase

        if verbose:
            print('channel_phase=', channel_phase)
            print('aos=', aos)
            print('lc=', lcs)
            print('amp=', amps)
            print('trend=', trends)

        # fit constant part of lc and trend
        #x, lcs_fit, trends_fit, readouts, sigmas = fit_eclipse_orbit(normal_orbit, aos, lcs, amps, channel_phase, shortFlag=use_short_states, longFlag=use_long_states)
        #trends_fit = numpy.zeros(n_pix) # just to illustrate difference

        # directly compute constant part of lc and trend for averaged eclipse data points
        trends_lin, lcs_lin = compute_trend(normal_orbit, aos, amps, channel_phase, **kwargs)

        #
        # write to db
        #

        vd = VarDark(normal_orbit, aos, lcs_lin, amps, trends_lin, channel_phase)
        self.fill_mtbl(vd)
        if self.created:
            self.append(vd)
        else:
            self.create(vd)

#- main ------------------------------------------------------------------------

use_short_states = False
use_long_states = True
vdd = VarDarkDb(verbose=False) # args=None
for orbit in range(27085,27195):
    vdd.calc_and_store_orbit(orbit, verbose=False, shortFlag=use_short_states, longFlag=use_long_states)
