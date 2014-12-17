#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module containing functionality for computation and storage of vardark (variable dark current of SCIA channel 8). 
"""

from __future__ import print_function, division

import h5py
import numpy as np
from numpy import cos, pi
from kapteyn import kmpfit
import logging

from ranges import remove_overlap, is_in_range, merge_ranges
from envisat import PhaseConverter
from sciamachy_module import petcorr, orbitfilter, get_darkstateid, read_extracted_states_ch8, NonlinCorrector
from scia_dark_functions import scia_dark_fun1, scia_dark_fun2

#- globals ---------------------------------------------------------------------

# orbit phase at which trending point lies. eclipse phase definition
trending_phase = 0.12

nlc = NonlinCorrector()
fname = '/SCIA/SDMF31/sdmf_extract_calib.h5'
n_pix = 1024

#- functions -------------------------------------------------------------------

def scia_dark_residuals1(p, data):
    x, y, yerr = data 
    return y - scia_dark_fun1(p, x)

def scia_dark_residuals1e(p, data):
    x, y, yerr = data 
    return (y - scia_dark_fun1(p, x)) / yerr

# should be raised if there is a truly unexpected error in the dark fit (data looks ok, but still failed)
class FitFailedError(Exception):
    pass

# should be raised if there is a truly unexpected error in the dark fit (data looks)
class NotEnoughDataError(Exception):
    pass

def fit_monthly(alldarks, orbit, verbose=False, kappasigma=False, debug_pixnr=None, short=False, give_errors=False, **kwargs):
    """
    fit dark model to two neighbouring monthly calibration orbits

    Parameters
    ----------

    alldarks : AllDarks object
        used to get underlying dark states from sdmf database in an easy and unambiguous way
    orbit : int
        absolute orbit number
    verbose : bool, optional
        if True, provide verbose output
    give_errors: bool, optional
        if True, return errors of fit parameters
    kappasigma : bool, optional
        if True, use kappasigma filter 
    """

    # minimum nr of degrees of freedom
    min_degrees = 8

    orbit_range = orbit-2., orbit+2.
    if verbose:
        print(orbit_range)

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = alldarks.get_range(orbit_range)

    #
    # fit it
    #

    x = ephases-orbit, pet #, coadd
    if verbose:
        print(pet.shape)
        print(pet)
        print(coadd)
    
    # small pets should not be taken veryseriously, unless specified or the rest is saturated..
    if not short:
        idx = pet < .49
        if np.sum(idx) > 0:
            idx = np.where(idx)[0]
            all_sigmas[idx,:] *= 50

    aoinfo = dict(fixed=False, limits=[0,10000])
    dcinfo = dict(fixed=False, limits=[-10000.,+500000.])
    amp1info = dict(fixed=False, limits=[-1000,+1000])
    trendinfo = dict(fixed=False, limits=[-1000,+1000])
    amp2info = dict(fixed=False, limits=[-1.,+1.])
    phase1info = dict(fixed=False, limits=[-3.,+3.])
    phase2info = dict(fixed=False, limits=[-3.,+3.])
    parinfo = [aoinfo,dcinfo,amp1info,trendinfo,phase1info,amp2info,phase2info]

    # prepare initial parameters
    ao0 = 3000
    dc0 = 4000
    amp1_0 = 10
    trend0 = 0
    amp2_0 = 0.1
    phase_offset1 = 0
    phase_offset2 = 0
    p0 = np.array([ao0, dc0, amp1_0, trend0, phase_offset1, amp2_0, phase_offset2])

    # ..and fit
    n_done = 0
    statuses = np.zeros(n_pix)
    res_phases = np.zeros(n_pix)
    res_phases2 = np.zeros(n_pix)
    for pixnr in range(n_pix):
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        idx_fin = np.isfinite(pix_readouts)
        if np.sum(idx_fin) > min_degrees:
            pix_readouts = pix_readouts[idx_fin]
            pix_sigmas = pix_sigmas[idx_fin]
            x = (ephases-orbit)[idx_fin], pet[idx_fin]

            # pass a
            idx = pix_sigmas == 0
            if np.sum(idx) > 0:
                pix_sigmas[idx] = 1000 # prohibit 0 sigmas as these crash kmpfit, just make it a huge sigma.
            fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)
            residual = scia_dark_fun2(fitobj.params, x) - pix_readouts
            dev_residual = np.std(residual)
            avg_residual = np.mean(residual)

            if kappasigma:
                # pass b : coarse kappa sigma filter 
                idx = (residual-avg_residual)**2 > 10.*dev_residual
                if np.sum(idx) > 0:
                    pix_sigmas[idx] *= 50
                fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)
                residual = scia_dark_fun2(fitobj.params, x) - pix_readouts
                dev_residual = np.std(residual)
                avg_residual = np.mean(residual)

                # pass c : finer kappa sigma filter 
                idx = (residual-avg_residual)**2 > 5.*dev_residual
                if np.sum(idx) > 0:
                    pix_sigmas[idx] *= 50
                fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)

            n_done += 1
        else:
            continue
        if verbose:
            print(pixnr, dev_residual, fitobj.params[4] % 1, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_phases[pixnr] = fitobj.params[4]
        res_phases2[pixnr] = fitobj.params[6]
        if (fitobj.status <= 0):
           raise FitFailedError('Error message = '+fitobj.message)

    # center the negative phase shift
    idx_neg = res_phases > 0.5
    res_phases[idx_neg] -= 0.5
    idx_neg = res_phases2 > 0.5
    res_phases2[idx_neg] -= 0.5
    # compute channel phase shift
    channel_phase1 = np.median(res_phases) % 1
    channel_phase2 = np.median(res_phases2) % 1
    np.sort(res_phases)

    if verbose:
        print('channel median phase =', channel_phase1, channel_phase2)
    phase1info = dict(fixed=True, limits=[-3.,+3.])
    phase2info = dict(fixed=True, limits=[-3.,+3.])
    parinfo = [aoinfo,dcinfo,amp1info,trendinfo,phase1info,amp2info,phase2info]  
    p0[4] = channel_phase1
    p0[6] = channel_phase2

    #
    # pass 2 - fix phase shift
    #

    aos = np.zeros(n_pix)
    lcs = np.zeros(n_pix)
    amps = np.zeros(n_pix)
    amps2 = np.zeros(n_pix)
    trends = np.zeros(n_pix)
    err_aos = np.zeros(n_pix)
    err_lcs = np.zeros(n_pix)
    err_amps = np.zeros(n_pix)
    err_amps2 = np.zeros(n_pix)
    err_phase1 = np.zeros(n_pix)
    err_trends = np.zeros(n_pix)
    err_phase2 = np.zeros(n_pix)
    statuses = np.zeros(n_pix)
    for pixnr in range(n_pix):
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        idx_fin = np.isfinite(pix_readouts)
        if np.sum(idx_fin) > min_degrees:
            pix_readouts = pix_readouts[idx_fin]
            pix_sigmas = pix_sigmas[idx_fin]
            x = (ephases-orbit)[idx_fin], pet[idx_fin]

            # pass a
            idx = pix_sigmas == 0
            if np.sum(idx) > 0:
                pix_sigmas[idx] = 1000 # prohibit 0 sigmas as these crash kmpfit, just make it a huge sigma.
            fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)
            residual = scia_dark_fun2(fitobj.params, x) - pix_readouts
            dev_residual = np.std(residual)
            avg_residual = np.mean(residual)

            if kappasigma:
                # pass b : coarse kappa sigma filter 
                idx = (residual-avg_residual)**2 > 10.*dev_residual
                if np.sum(idx) > 0:
                    pix_sigmas[idx] *= 50
                fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)
                residual = scia_dark_fun2(fitobj.params, x) - pix_readouts
                dev_residual = np.std(residual)
                avg_residual = np.mean(residual)

                # pass c : finer kappa sigma filter 
                idx = (residual-avg_residual)**2 > 5.*dev_residual
                if np.sum(idx) > 0:
                    pix_sigmas[idx] *= 50
                fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)

            n_done += 1
        else:
            continue
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
        err_aos[pixnr] = fitobj.xerror[0]
        err_lcs[pixnr] = fitobj.xerror[1]
        err_amps[pixnr] = fitobj.xerror[2]
        err_trends[pixnr] = fitobj.xerror[3]
        err_phase1[pixnr] = fitobj.xerror[4]
        err_amps2[pixnr] = fitobj.xerror[5]
        err_phase2[pixnr] = fitobj.xerror[6]

    channel_amp2 = np.median(amps2[np.where(statuses > 0)])

    if debug_pixnr is not None:
        import matplotlib.pyplot as plt
        plt.cla()
        p = aos[debug_pixnr], lcs[debug_pixnr], amps[debug_pixnr], trends[debug_pixnr], channel_phase1, amps2[debug_pixnr], channel_phase2
        x = ephases-orbit, pet
        model = scia_dark_fun2(p,x)
        print("model=", model)
        plt.plot(ephases-orbit, all_readouts[:,debug_pixnr], 'bo', label="data")
        plt.plot(ephases-orbit, model, 'go', label="model")
        plt.legend(loc="best")
        plt.show()

    if give_errors:
        errors = {"aos":err_aos, 
                  "off":err_lcs, 
                  "amps":err_amps, 
                  "amps2":np.median(err_amps2), 
                  "phase1":np.median(err_phase1), 
                  "phase2":np.median(err_phase2), 
                  "trends":err_trends}
        return channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends, errors
    else:
        return channel_phase1, channel_phase2, aos, lcs, amps, channel_amp2, trends

def fit_eclipse_orbit(alldarks, orbit, aos, lcs, amps, amp2, channel_phaseshift, channel_phaseshift2, 
                      give_errors=False, verbose=False, short=False, kappasigma=False, **kwargs):
    """
    fit dark model to an orbits with normal (eclipse) dark states
    use parameters computed from nearest calibration orbits

    Parameters
    ----------

    alldarks : AllDarks object
        used to get underlying dark states from sdmf database in an easy and unambiguous way
    orbit : int
        absolute orbit number
    aos : numpy array, 1024 floats 
        analog offset per pixel for this orbit
    lcs : numpy array, 1024 floats 
        leakage current (constant part) per pixel for this orbit
    amps : numpy array, 1024 floats 
        orbital variation amplitude per pixel for this orbit
    channel_amp2 : float
        amplitude of first harmonic of orbital variation (channel mean)
    channel_phaseshift : float
        phase shift of orbital variation (channel mean)
    channel_phaseshift2 : float
        phase shift of first harmonic orbital variation (channel mean)
    give_errors : bool, optional
        if set, return also fit uncertainty
    verbose : bool, optional
        if set, provide verbose output
    kappasigma : bool, optional
        if set, use kappasigma filter 
    """

    # minimum nr of degrees of freedom
    min_degrees = 3

    #
    # get all dark data
    # 

    orbit_range = orbit-1.25, orbit+1.25
    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = alldarks.get_range(orbit_range)

    if ephases.size <= 2:
        raise NotEnoughDataError("Not enough datapoints for fit: "+str(ephases.size))

    #
    # initialize data points (coadd is always 1 as we only do ch8 nadir)
    #

    x = ephases-orbit, pet #, coadd

    # small pets should not be taken veryseriously, unless specified or the rest is saturated..
    if not short:
        idx = pet < .49
        if np.sum(idx) > 0:
            idx = np.where(idx)[0]
            all_sigmas[idx,:] *= 50

    #
    # setup attributes of fit parameters
    #

    # note the limits are just slightly wider than in the monthly fit. we do this to get rid of float32->float64 conversion errors!
    aoinfo = dict(fixed=True, limits=[-0.1,11000.1])
    dcinfo = dict(fixed=False, limits=[-10000.1,+500000.1])
    amp1info = dict(fixed=True, limits=[-1000.1,+1000.1])
    trendinfo = dict(fixed=False, limits=[-1000.1,+1000.1])
    amp2info = dict(fixed=True, limits=[-1.01,+1.01])
    phase1info = dict(fixed=True, limits=[-3.01,+3.01])
    phase2info = dict(fixed=True, limits=[-3.01,+3.01])
    parinfo = [aoinfo,dcinfo,amp1info,trendinfo,phase1info,amp2info,phase2info] 

    #
    # setup result data arrays
    #

    n_done = 0
    statuses = np.zeros(n_pix)
    res_trends = np.empty(n_pix)
    res_lcs = np.empty(n_pix)
    res_trends[:] = np.nan
    res_lcs[:] = np.nan
    err_trends = np.empty(n_pix)
    err_lcs = np.empty(n_pix)
    err_trends[:] = np.nan
    err_lcs[:] = np.nan
    uncertainty = np.empty(n_pix)
    uncertainty[:] = np.nan

    #
    # ..and fit
    #

    for pixnr in range(n_pix):
        # prepare initial parameters.. 
        p0 = np.array([aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phaseshift, amp2, channel_phaseshift2]) 
        pix_readouts = all_readouts[:,pixnr]
        pix_sigmas = all_sigmas[:,pixnr]
        idx_fin = np.isfinite(pix_readouts)
        if np.sum(idx_fin) > min_degrees:
            pix_readouts = pix_readouts[idx_fin]
            pix_sigmas = pix_sigmas[idx_fin]
            x = (ephases-orbit)[idx_fin], pet[idx_fin]

            #print(orbit, pixnr, p0, parinfo)
            idx = pix_sigmas == 0
            if np.sum(idx) > 0:
                pix_sigmas[idx] = 1000 # prohibit 0 sigmas as these crash kmpfit, just make it a huge sigma.
            fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, xtol=1e-8, parinfo=parinfo)
            residual = scia_dark_fun2(fitobj.params, x) - pix_readouts
            dev_residual = np.std(residual)
            avg_residual = np.mean(residual)

            if kappasigma:
                # pass b : coarse kappa sigma filter 
                idx = (residual-avg_residual)**2 > 10.*dev_residual
                if np.sum(idx) > 0:
                    pix_sigmas[idx] *= 50
                fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)
                residual = scia_dark_fun2(fitobj.params, x) - pix_readouts
                dev_residual = np.std(residual)
                avg_residual = np.mean(residual)

                # pass c : finer kappa sigma filter 
                idx = (residual-avg_residual)**2 > 5.*dev_residual
                if np.sum(idx) > 0:
                    pix_sigmas[idx] *= 50
                fitobj = kmpfit.simplefit(scia_dark_fun2, p0, x, pix_readouts, err=pix_sigmas, ftol=1e-8, parinfo=parinfo)

            n_done += 1
        else:
            continue
        if verbose:
            print(pixnr, fitobj.message)
            print(p0, x, pix_readouts, pix_sigmas)
        if (fitobj.status <= 0):
            raise FitFailedError(fitobj.message)
        statuses[pixnr] = fitobj.status
        res_lcs[pixnr] = fitobj.params[1]
        err_lcs[pixnr] = fitobj.stderr[1]
        res_trends[pixnr] = fitobj.params[3]
        err_trends[pixnr] = fitobj.stderr[3]
        uncertainty[pixnr] = np.std(scia_dark_fun2(fitobj.params, x) - pix_readouts)

    if give_errors:
        return x, res_lcs, res_trends, err_lcs, err_trends, all_readouts, all_sigmas, uncertainty
    else:
        return x, res_lcs, res_trends, all_readouts, all_sigmas

def read_ch8_darks(orbit_range, stateid):
    """
    Read reduced channel 8dark states with specified state id from sdmf_extract_calib.
    Corrects non-linearity.

    Parameters
    ----------
    orbit_range: 2 element list of ints 
        list comprising start and end orbit
    stateid: int
        state id of the dark

    Returns
    -------
    jds : arraylike, float64
        julian days of state executions
    readouts : arraylike, 2d, float
        readouts
    noise : arraylike, 2d, float
        measure of uncertainty of the readouts
    tdet : float
        detector temperature
    pet : float
        pixel exposure time (in seconds)
    coadd : int
        co-adding factor
    """
    states = read_extracted_states_ch8(orbit_range, stateid, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states['mtbl']
    jds = state_mtbl['julianDay'][:]
    readouts = states['readoutMean']
    noise = states['readoutNoise']
    n_exec = readouts.shape[0]
    for idx_exec in range(n_exec):
        readouts[idx_exec,:] = nlc.correct_ch8(readouts[idx_exec,:])
    pet = np.zeros(n_exec) + states['pet'][0]
    coadd = np.zeros(n_exec) + states['coadd'][0]
    noise = noise / np.sqrt(coadd[0])
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
        self.jds_ = np.array([])
        self.jds_ = np.array([])
        self.readouts_ = np.empty([0, n_pix])
        self.noise_ = np.empty([0, n_pix])
        self.tdet_ = np.array([])
        self.pet_ = np.array([])
        self.coadd_ = np.array([])
        self.ephases = np.array([])
        self.range_list = []
        self.pc = PhaseConverter()
        return

    def _register_range(self, orbit_range):
        """
        register orbit range and automatically remove overlap
        """
        self.range_list.append((orbit_range[0], orbit_range[1])) # make sure we're entering tuples into the list!
        self.range_list = remove_overlap(self.range_list)
        self.range_list = merge_ranges(self.range_list)
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
            if ((orbit_range[0] < 43361) and (orbit_range[1] > 43362)):
                raise Exception("orbit range "+str(orbit_range)+" crosses OCR43!") # dark state definition should not change during 1 orbit range
            orbit = (orbit_range[0]+orbit_range[1]) // 2 # take the center to be sure we're not on the other side of OCR43
            # a bit hacky, but since we really do not have these orbits with 1.0s PET, it's best to just skip instead of some complicated exception handling 
            if ((petje == 1.0) or (petje == 0.125)) and (orbit < 4151):
                continue
            stateid = get_darkstateid(petje, orbit)
            print(petje, stateid)
            jds_, readouts_, noise_, tdet_, pet_, coadd_ = read_ch8_darks(orbit_range, stateid)
            self.jds_ = np.concatenate((self.jds_, jds_))
            self.readouts_ = np.concatenate((self.readouts_, readouts_))
            self.noise_ = np.concatenate((self.noise_, noise_))
            self.tdet_ = np.concatenate((self.tdet_, tdet_))
            self.pet_ = np.concatenate((self.pet_, petje+np.zeros(jds_.size)))
            self.coadd_ = np.concatenate((self.coadd_, coadd_))

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
        self.jds_, idx = np.unique(self.jds_, return_index=True)
        self.readouts_ = self.readouts_[idx,:]
        self.noise_ = self.noise_[idx,:]
        self.tdet_ = self.tdet_[idx]
        self.pet_ = self.pet_[idx]
        self.coadd_ = self.coadd_[idx]

        # get eclipse phases + orbits
        ephases, orbits = self.pc.get_phase(self.jds_, getOrbits=True)
        ephases += orbits

        # filter out sunrise-affected data
        ephases1 = np.mod(ephases, 1.)
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
                print("lump from 43362")
                self._lumpr([43362, last_orbit])
            else:
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
        return np.sum(idx), 0, self.pet[idx], self.coadd[idx], self.readouts[idx,:], self.noise[idx,:], self.ephases[idx] 

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

def load_sdmf30_dark(orbit):
    """
    load channel 8 analog offset and dark current from SDMF3.0 product: a good fall-back for nearly saturated pixels.
    """
    fname = "/SCIA/SDMF30/sdmf_dark.h5"
    fid = h5py.File(fname, "r")
    orbits = fid["orbitList"]
    idx = orbits[:] == orbit
    if np.sum(idx) == 0:
        raise Exception("orbit ["+str(orbit)+"] not found in ["+fname+"]")
    idx = np.where(idx)[0][0]
    dset_ao = fid["analogOffset"]
    dset_dc = fid["darkCurrent"]
    aos = (dset_ao[7*1024:,idx]).flatten()
    dcs = (dset_dc[7*1024:,idx]).flatten()
    return aos, dcs

#- main ------------------------------------------------------------------------

if __name__ == "__main__":
    """
    just a test..
    """

    of = orbitfilter()
    orbit = of.get_closest_monthly(43500)
    print(orbit)

    ad = AllDarks([0.125, 0.5, 1.0])

    ret = fit_monthly(ad, orbit, verbose=False, kappasigma=False, debug_pixnr=489, short=False, give_errors=True)
    channel_phase1, channel_phase2, aos, dcs, amps, channel_amp2, trends = ret 

    pixels = [399,406,415,423,424,426,431,433,458,463,465,484,504,514,517,532,534,544,549,562,563,574,577,578,582,598,600,604,613,615,597]

    print("channel_phase1=",channel_phase1)
    print("channel_phase2=",channel_phase2)
    print("channel_amp2=",channel_amp2)
    for pixnr in pixels:
        print("PIXNR", pixnr)
        print("ao=",aos[pixnr])
        print("dc=",dcs[pixnr])
        print("amp=",amps[pixnr])
        print("trend=",trends[pixnr])
