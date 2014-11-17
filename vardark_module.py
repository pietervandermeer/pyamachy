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

from __future__ import print_function, division

import h5py
import numpy
import numpy as np
from numpy import cos, pi
from kapteyn import kmpfit
import matplotlib.pyplot as plt

from ranges import remove_overlap, is_in_range, merge_ranges
from envisat import PhaseConverter
from sciamachy_module import NonlinCorrector, read_extracted_states, petcorr, orbitfilter, get_darkstateid, read_extracted_states_
from scia_dark_functions import scia_dark_fun1, scia_dark_fun2

#-------------------------SECTION VERSION-----------------------------------

_swVersion = {'major': 0,
              'minor': 3,
              'revision' : 0}
_calibVersion = {'major': 0,
                 'minor': 1,
                 'revision' : 0}
_dbVersion = {'major': 0,
              'minor': 2,
              'revision' : 0}

#- globals ---------------------------------------------------------------------

# orbit phase at which trending point lies. eclipse phase definition
trending_phase = 0.12 # orb < 43362
#trending_phase = 0.35

nlc = NonlinCorrector()
phaseconv = PhaseConverter()
fname = '/SCIA/SDMF31/sdmf_extract_calib.h5'
n_pix = 1024

#- functions -------------------------------------------------------------------

def scia_dark_residuals1(p, data):
    x, y, yerr = data 
    return y - scia_dark_fun1(p, x)

def scia_dark_residuals1e(p, data):
    x, y, yerr = data 
    return (y - scia_dark_fun1(p, x)) / yerr

# fit dark model to two neighbouring monthly calibration orbits
#def fit_monthly(orbit, shortFlag=False, longFlag=True):
def fit_monthly(alldarks, orbit, verbose=False, **kwargs):

    orbit_range = orbit-.5, orbit+2.5
    if verbose:
        print(orbit_range)

    #n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = extract_dark_states(orbit, **kwargs)
    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = alldarks.get_range(orbit_range)

    #
    # fit it
    #

    x = ephases-orbit, pet #, coadd
    if verbose:
        print(pet.shape)
        print(pet)
        print(coadd)
    
    aoinfo = dict(fixed=False, limits=[0,10000])
    lcinfo = dict(fixed=False, limits=[-100000.,+100000.])
    amp1info = dict(fixed=False, limits=[-1000,+1000])
    trendinfo = dict(fixed=False, limits=[-1000,+1000])
    amp2info = dict(fixed=False, limits=[-1.,+1.])
    phase1info = dict(fixed=False, limits=[-3.,+3.])
    phase2info = dict(fixed=False, limits=[-3.,+3.])
#    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info] 
    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info,amp2info,phase2info]

    # prepare initial parameters
    ao0 = 4000
    dc0 = 4000
    amp1_0 = 10
    trend0 = 0
    amp2_0 = 0.1
    phase_offset1 = 0
    phase_offset2 = 0
    #p0 = numpy.array([ao0, dc0, amp1_0, trend0, phase_offset1])
    p0 = numpy.array([ao0, dc0, amp1_0, trend0, phase_offset1, amp2_0, phase_offset2])

    # ..and fit
    n_done = 0
    statuses = numpy.zeros(n_pix)
    res_phases = numpy.zeros(n_pix)
    res_phases2 = numpy.zeros(n_pix)
    for pixnr in range(n_pix):
        #print(pixnr)
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        #print('readouts = ', pix_readouts)
        #print('sigmas = ', pix_sigmas)
        #print(numpy.sum(pix_sigmas))
        if not numpy.isnan(numpy.sum(pix_readouts)) and numpy.all(pix_sigmas != 0):
#            fitobj = kmpfit.simplefit(scia_dark_fun1, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            n_done += 1
        else:
            continue
        if verbose:
            print(pixnr, fitobj.params[4] % 1, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_phases[pixnr] = fitobj.params[4]
        res_phases2[pixnr] = fitobj.params[6]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()

    # center the negative phase shift
    idx_neg = res_phases > 0.5
    res_phases[idx_neg] -= 0.5
    idx_neg = res_phases2 > 0.5
    res_phases2[idx_neg] -= 0.5
    # compute channel phase shift
    channel_phase1 = numpy.median(res_phases[numpy.where(statuses > 0)] ) % 1
    channel_phase2 = numpy.median(res_phases2[numpy.where(statuses > 0)] ) % 1
    if verbose:
        print('channel median phase =', channel_phase1, channel_phase2)
    phase1info = dict(fixed=True, limits=[-3.,+3.])
    phase2info = dict(fixed=True, limits=[-3.,+3.])
    #parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info]
    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info,amp2info,phase2info]  
    p0[4] = channel_phase1
    p0[6] = channel_phase2

    #
    # pass 2 - fix phase shift
    #

    aos = numpy.zeros(n_pix)
    lcs = numpy.zeros(n_pix)
    amps = numpy.zeros(n_pix)
    amps2 = numpy.zeros(n_pix)
    trends = numpy.zeros(n_pix)
    statuses = numpy.zeros(n_pix)
    for pixnr in range(n_pix):
        #print(pixnr)
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        #print('readouts = ', pix_readouts)
        #print('sigmas = ', pix_sigmas)
        #print(numpy.sum(pix_sigmas))
        if not numpy.isnan(numpy.sum(pix_readouts)) and numpy.all(pix_sigmas != 0):
            #fitobj = kmpfit.simplefit(scia_dark_fun1, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
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
        amps2[pixnr] = fitobj.params[5]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()

    channel_amp2 = numpy.median(amps2[numpy.where(statuses > 0)])

    # TODO: initial outlier removal pass?

    #plt.cla()
    #plt.scatter(all_state_phases, all_readouts[:,pixnr])
    #plt.show()

    # TODO: return refinement parameters as well.

    return channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends

def fit_eclipse_orbit(alldarks, orbit, aos, lcs, amps, amp2, channel_phaseshift, channel_phaseshift2, 
                      give_errors=False, verbose=False, **kwargs):
    """
    fit dark model to an orbits with normal (eclipse) dark states
    use parameters computed from nearest calibration orbits
    """

    orbit_range = orbit, orbit+2
    #n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = extract_dark_states(orbit, **kwargs)
    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = alldarks.get_range(orbit_range)

    if ephases.size <= 2:
        raise Exception("Not enough datapoints for fit:", ephases.size)

    #
    # fit it
    #

    x = ephases-orbit, pet #, coadd

    # note the limits are just slightly wider than in the monthly fit. we do this to get rid of float32->float64 conversion errors!
    aoinfo = dict(fixed=True, limits=[-0.1,10000.1])
    lcinfo = dict(fixed=False, limits=[-100000.1,+100000.1])
    amp1info = dict(fixed=True, limits=[-1000.1,+1000.1])
    trendinfo = dict(fixed=False, limits=[-1000.1,+1000.1])
    amp2info = dict(fixed=True, limits=[-1.01,+1.01])
    phase1info = dict(fixed=True, limits=[-3.01,+3.01])
    phase2info = dict(fixed=True, limits=[-3.01,+3.01])
    parinfo = [aoinfo,lcinfo,amp1info,trendinfo,phase1info,amp2info,phase2info] 

    # ..and fit
    n_done = 0
    statuses = numpy.zeros(n_pix)
    res_trends = numpy.empty(n_pix)
    res_lcs = numpy.empty(n_pix)
    res_trends[:] = numpy.nan
    res_lcs[:] = numpy.nan
    err_trends = numpy.empty(n_pix)
    err_lcs = numpy.empty(n_pix)
    err_trends[:] = numpy.nan
    err_lcs[:] = numpy.nan
    uncertainty = numpy.empty(n_pix)
    uncertainty[:] = numpy.nan
    for pixnr in range(n_pix):
        # prepare initial parameters.. 
        p0 = numpy.array([aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phaseshift, amp2, channel_phaseshift2]) 
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        if (aos[pixnr] is not 0) and (not numpy.isnan(numpy.sum(pix_readouts))) and (numpy.all(pix_sigmas != 0)) and (x[0].size > 0):
            #print(orbit, pixnr, x)
            #print(orbit, pixnr, p0, parinfo)
            fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            n_done += 1
        else:
            continue
        if verbose:
            print(pixnr, fitobj.message)
            print(p0, x, pix_readouts, pix_sigmas)
        if (fitobj.status <= 0):
            raise Exception(fitobj.message)
        statuses[pixnr] = fitobj.status
        res_lcs[pixnr] = fitobj.params[1]
        err_lcs[pixnr] = fitobj.stderr[1]
        res_trends[pixnr] = fitobj.params[3]
        err_trends[pixnr] = fitobj.stderr[3]
        uncertainty[pixnr] = np.std(scia_dark_fun2(fitobj.params, x) - pix_readouts)
        #else:
        #   print("Optimal parameters: ", fitobj.params)

    if give_errors:
        return x, res_lcs, res_trends, err_lcs, err_trends, all_readouts, all_sigmas, uncertainty
    else:
        return x, res_lcs, res_trends, all_readouts, all_sigmas

def read_ch8_darks(orbit_range, stateid):
    states = read_extracted_states_(orbit_range, stateid, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states['mtbl']
    jds = state_mtbl['julianDay'][:]
    readouts = states['readoutMean']
    noise = states['readoutNoise']
    n_exec = readouts.shape[0]
    for idx_exec in range(n_exec):
        readouts[idx_exec,:] = nlc.correct_ch8(readouts[idx_exec,:])
    pet = numpy.zeros(n_exec) + states['pet'][0]
    coadd = numpy.zeros(n_exec) + states['coadd'][0]
    noise = noise / numpy.sqrt(coadd[0])
    phases = state_mtbl['orbitPhase'][:]
    tdet = state_mtbl['detectorTemp'][:]
    tdet = tdet[:,7].flatten() # ch8
    return jds, readouts, noise, tdet, pet, coadd

class AllDarks():
    """
    Class that buffers ranges of darks for the Pixel Exposure Times (PETs) it was initialized for,
    and returns dark data for specified range all from buffer.
    It is relatively smart in the sense that it handles OCR43 and bad orbit slicing automatically.
    It also requires minimal initialization. 
    """
    def __init__(self, petlist):
        """
        Parameters
        ----------

        petlist : list or tuple of floats
            PETs of the dark states
        """
        self.petlist = petlist
        self.jds_ = numpy.array([])
        self.jds_ = numpy.array([])
        self.readouts_ = numpy.empty([0, n_pix])
        self.noise_ = numpy.empty([0, n_pix])
        self.tdet_ = numpy.array([])
        self.pet_ = numpy.array([])
        self.coadd_ = numpy.array([])
        self.ephases = numpy.array([])
        self.range_list = []

    def _register_range(self, orbit_range):
        """
        register orbit range and automatically remove overlap
        """
        self.range_list.append((orbit_range[0], orbit_range[1])) # make sure we're entering tuples into the list!
        print("post append:", self.range_list)
        self.range_list = remove_overlap(self.range_list)
        print("post remove_overlap:", self.range_list)
        self.range_list = merge_ranges(self.range_list)
        print("post merge:", self.range_list)
        return

    def is_registered(self, orbit_range):
        """
        check if the range is within the previously fetched ranges.
        """
        return is_in_range(self.range_list, orbit_range)

    def _lump(self, orbit_range):
        """
        a simple approach at fetching data. just append everything and let finalize() deduplicate.
        """
        # concat all states
        for petje in self.petlist:
            #print("PET=",petje)
            stateid = get_darkstateid(petje, orbit_range[0]) # TODO: this could be on the border of two state definitions
            jds_, readouts_, noise_, tdet_, pet_, coadd_ = read_ch8_darks(orbit_range, stateid)
            #print(stateid, pet_)
            self.jds_ = numpy.concatenate((self.jds_, jds_))
            self.readouts_ = numpy.concatenate((self.readouts_, readouts_))
            self.noise_ = numpy.concatenate((self.noise_, noise_))
            self.tdet_ = numpy.concatenate((self.tdet_, tdet_))
            self.pet_ = numpy.concatenate((self.pet_, pet_))
            self.coadd_ = numpy.concatenate((self.coadd_, coadd_))

        return

    def _lumpl(self, orbit_range):
        """
        lump with left side extended by one orbit. (because of orbit slicing in sdmf extract)
        """
        self._register_range(orbit_range)
        self._lump([orbit_range[0]-1, orbit_range[1]])

    def _lumpr(self, orbit_range):
        """
        lump with right side extended by one orbit. (because of orbit slicing in sdmf extract)
        """
        self._register_range(orbit_range)
        self._lump([orbit_range[0], orbit_range[1]+1])

    def _lumplr(self, orbit_range):
        """
        lump with left and right side extended by one orbit. (because of orbit slicing in sdmf extract)
        """
        self._register_range(orbit_range)
        self._lump([orbit_range[0]-1, orbit_range[1]+1])

    def _finalize(self):
        """
        finalizes (deduplicates) self.jds_, self.readouts_ .. etc to self.jds, self.readouts, etc
        eclipse phases are computed and sunrise is filtered out
        """

        # deduplicate after having lumped
        self.jds_, idx = numpy.unique(self.jds_, return_index=True)
        self.readouts_ = self.readouts_[idx,:]
        self.noise_ = self.noise_[idx,:]
        self.tdet_ = self.tdet_[idx]
        self.pet_ = self.pet_[idx]
        self.coadd_ = self.coadd_[idx]

        # get eclipse phases + orbits
        ephases, orbits = phaseconv.get_phase(self.jds_, getOrbits=True)
        ephases += orbits

        # filter out sunrise-affected data
        ephases1 = numpy.mod(ephases, 1.)
        idx_nosunrise = (ephases1 < .35) | (ephases1 > .42)
        self.jds = self.jds_[idx_nosunrise]
        self.ephases = ephases[idx_nosunrise]
        self.readouts = self.readouts_[idx_nosunrise,:]
        self.noise = self.noise_[idx_nosunrise,:]
        self.tdet = self.tdet_[idx_nosunrise]
        self.pet = self.pet_[idx_nosunrise]
        self.coadd = self.coadd_[idx_nosunrise]

        return

    def buffer_range(self, orbit_range):
        """
        buffers the specified orbit range of darks.
        takes into account OCR43 dark state definition change. 
        takes into account orbit slicing issue in level 0 and sdmf extract db

        Parameters
        ----------

        orbit_range : tuple or list, 2 integerss
            orbit range
        """
        if not self.is_registered(orbit_range):
            first_orbit = orbit_range[0]
            last_orbit = orbit_range[1]
            if first_orbit < 43362 and last_orbit >= 43362:
                print("lump upto 43361")
                self._lumpl([first_orbit, 43361])
                print("LO???",self.range_list)
                print("lump from 43362")
                self._lumpr([43362, last_orbit])
                print("HI???",self.range_list)
            else:
                print("buffer_range(): lumplr", orbit_range)
                self._lumplr(orbit_range)
            self._finalize()

    def get_range(self, orbit_range, autoLump=True):
        """
        retrieves the specified orbit range of darks for the caller. 
        takes into account OCR43 dark state definition change. 
        takes into account orbit slicing issue in level 0 and sdmf extract db

        Parameters
        ----------

        orbit_range : tuple or list, 2 integerss
            orbit range
        autoLump : boolean, optional
            set this to automatically buffer darks that are not yet in the buffer
        """
        if autoLump:
            self.buffer_range(orbit_range)
        idx = (self.ephases.astype('i') >= orbit_range[0]) & (self.ephases.astype('i') <= orbit_range[1])
        return numpy.sum(idx), 0, self.pet[idx], self.coadd[idx], self.readouts[idx,:], self.noise[idx,:], self.ephases[idx] 

# compute thermal background trend between two normal orbits 
# also computes actual thermal background offset (excluding trend or oscillation)
def compute_trend(alldarks, orbit, aos, amps, amp2, phaseshift, phaseshift2, **kwargs):

    orbit_range = orbit-.5, orbit+2.5
    #print(orbit_range)

    trends = numpy.empty(n_pix)
    trends[:] = numpy.nan
    lcs = numpy.empty(n_pix)
    lcs[:] = numpy.nan

    #n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = extract_dark_states(orbit, **kwargs)
    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = alldarks.get_range(orbit_range)
    print(n_exec)
    if n_exec == 0:
        #return numpy.zeros(n_pix), numpy.zeros(n_pix)
        return trends, lcs

    #print(ephases)
    #print(pet)
    # divide by coadding factor, subtract analog offset and divide by exposure time to get thermal background in BU/s
    thermal_background = all_readouts
    #thermal_background /= numpy.matrix(coadd).T * numpy.ones(n_pix) # this was already done.. 
    thermal_background -= (numpy.matrix(aos).T * numpy.ones(n_exec)).T
    thermal_background /= numpy.matrix(pet).T * numpy.ones(n_pix)

    abs_dist1 = numpy.abs(ephases - (trending_phase+orbit))
    abs_dist2 = numpy.abs(ephases - (trending_phase+orbit+1))
    idx1 = (numpy.argsort(abs_dist1))[0:6]
    idx2 = (numpy.argsort(abs_dist2))[0:6]
    #print(idx2)
    phi1 = numpy.mean(ephases[idx1])
    phi2 = numpy.mean(ephases[idx2])
    #print(phi1,phi2)
    idx0 = numpy.argmin(ephases)
    phi0 = ephases[idx0]
    phi1_ = phi1 % 1

    # loop over all channel pixels
    for pixnr in range(n_pix):
        if aos[pixnr] is 0:
            continue
        if numpy.all(all_sigmas[:,pixnr] is 0):
            continue 
        background = thermal_background[:,pixnr]
        if numpy.isnan(numpy.sum(background)):
            continue 

        avg1 = numpy.mean(background[idx1])
        avg2 = numpy.mean(background[idx2])
        trends[pixnr] = (avg2-avg1)/(phi2-phi1)
        lcs[pixnr] = avg1 - trends[pixnr]*phi1_ - amps[pixnr]*( cos(2*pi*(phi1_+phaseshift)) + amp2*cos(4*pi*(phi1_+phaseshift2)) )
#        lcs[pixnr] = avg1 - trends[pixnr]*phi1_ - amps[pixnr]*( cos(2*pi*(phi1_+phaseshift)) + amp2*cos(4*pi*(phi1_+phaseshift2)) )

    return trends, lcs

# TODO: OUTDATED.. replaced with another class of the same name. up for removal
class VarDark:
    def __init__(self, orbit,ao, lc, amp, trend, phaseshift, amp2, phaseshift2):
        self.numPixels = 1024
        self.absOrbit = orbit
        self.ao = ao
        self.lc = lc
        self.amp = amp
        self.trend = trend
        self.phaseshift = phaseshift
        self.amp2 = amp2
        self.mtbl = dict()
        self.mtbl['julianDay'] = 1.234 # TODO convert orbit to JD
        self.mtbl['phaseShift2'] = phaseshift2
        self.mtbl['amp2'] = amp2
        self.mtbl['obmTemp'] = 999 # TODO
        self.mtbl['detTemp'] = numpy.array([999,999,999,999,999,999,999,999]) # TODO

class VarDarkDb:
    def __init__( self, args=None, db_name='./sdmf_vardark.h5',
                  truncate=False, calibration=None, verbose=False, allDarks=None ):
        self.monthly_orbit = -1
        self.created = False
        self.ofilt = orbitfilter()
        # initialize a minimal version
        self.alldarks = AllDarks([0.5, 1.0])

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
        return

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

        fmtMTBL  = 'float64,a20,uint16,uint16,float32,float32,float32,float32,8float32'
        nameMTBL = ('julianDay','entryDate','absOrbit','quality','phaseShift','phaseShift2','amp2','obmTemp','detTemp')

        self.mtbl = np.empty(1, dtype=fmtMTBL )
        self.mtbl.dtype.names = nameMTBL
        self.mtbl['julianDay'] = vardark.mtbl['julianDay']
        self.mtbl['entryDate'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.mtbl['absOrbit'] = vardark.absOrbit
        self.mtbl['quality'] = 0
        self.mtbl['phaseShift'] = vardark.phaseshift
        self.mtbl['obmTemp'] = vardark.mtbl['obmTemp']
        self.mtbl['detTemp'] = vardark.mtbl['detTemp']
        self.mtbl['phaseShift2'] = vardark.mtbl['phaseShift2']
        self.mtbl['amp2'] = vardark.mtbl['amp2']

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
        monthly_orbit = self.ofilt.get_closest_monthly(orbit)
        if self.monthly_orbit != monthly_orbit:
            # calculate if monthly if not buffered
            channel_phase, channel_phase2, aos, lcs, amps, amp2, trends = fit_monthly(self.alldarks, monthly_orbit, **kwargs)
            self.monthly_orbit = monthly_orbit
            self.aos = aos
            self.lcs = lcs
            self.amps = amps
            self.amp2 = amp2
            self.trends = trends
            self.channel_phase = channel_phase 
            self.channel_phase2 = channel_phase2
        else:
            # just use the buffered version
            aos = self.aos
            lcs = self.lcs
            amps = self.amps
            trends = self.trends
            channel_phase = self.channel_phase
            amp2 = self.amp2
            channel_phase2 = self.channel_phase2

        if verbose:
            print('channel_phase=', channel_phase)
            print('channel_phase2=', channel_phase2)
            print('aos=', aos)
            print('lc=', lcs)
            print('amp=', amps)
            print('amp2=', amp2)
            print('trend=', trends)

        # fit constant part of lc and trend
        #x, lcs_fit, trends_fit, readouts, sigmas = fit_eclipse_orbit(normal_orbit, aos, lcs, amps, channel_phase, shortFlag=use_short_states, longFlag=use_long_states)
        #trends_fit = numpy.zeros(n_pix) # just to illustrate difference

        # directly compute constant part of lc and trend for averaged eclipse data points
        trends_lin, lcs_lin = compute_trend(normal_orbit, aos, amps, amp2, channel_phase, channel_phase2, **kwargs)

        #
        # write to db
        #

        vd = VarDark(normal_orbit, aos, lcs_lin, amps, trends_lin, channel_phase, amp2, channel_phase2)
        self.fill_mtbl(vd)
        if self.created:
            self.append(vd)
        else:
            self.create(vd)

def load_varkdark_orbit(orbit, shortMode, give_uncertainty=False, fname=None):
    if fname is None:
        basename = "vardark"
        if shortMode:
            fname = basename+"_short.h5"
        else:
            fname = basename+"_long.h5"
    print(fname)
    fid = h5py.File(fname, "r")
    orbit_dset = fid["dim_orbit"]
    idx_fid = orbit_dset[:] == orbit
    if np.sum(idx_fid) == 0:
        raise Exception("no thermal background found for orbit "+str(orbit)+"!")
    therm_dset = fid["varDark"]
    if give_uncertainty:
        uncertainty_dset = fid["uncertainties"]
        uncertainty = uncertainty_dset[idx_fid,:]
    thermal_background = therm_dset[idx_fid,:,:]
    phases = fid["dim_phase"][:] + orbit
    fid.close()

    basename = "interpolated_monthlies"
    if shortMode:
        fname = basename+"_short.h5"
    else:
        fname = basename+"_long.h5"
    fin = h5py.File(fname, "r")
    in_orbitlist = fin["orbits"]
    idx_fin = in_orbitlist[:] == orbit
    if np.sum(idx_fin) == 0:
        raise Exception("no analog offset data found for orbit "+str(orbit)+"!")
    analog_offset = fin['aos'][idx_fin,:]
    fin.close()

    if give_uncertainty:
        return phases, orbit, thermal_background, analog_offset, uncertainty
    else:
        return phases, orbit, thermal_background, analog_offset

#- main ------------------------------------------------------------------------

# a test product generating entity. runs only when this module is run as a program.
if __name__ == "__main__":
    start_orbit = 27085
    end_orbit = 27195
    print("building vardark long..")
    vddl = VarDarkDb(verbose=False, db_name="sdmf_vardark_long.h5") # args=None
    for orbit in range(start_orbit, end_orbit):
        vddl.calc_and_store_orbit(orbit, verbose=False, shortFlag=False, longFlag=True)
    print("building vardark short..")
    vdds = VarDarkDb(verbose=False, db_name="sdmf_vardark_short.h5") # args=None
    for orbit in range(start_orbit, end_orbit):
        vdds.calc_and_store_orbit(orbit, verbose=False, shortFlag=True, longFlag=False)

