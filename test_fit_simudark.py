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

import numpy
from numpy import cos, pi
from kapteyn import kmpfit
import matplotlib.pyplot as plt

from sciamachy_module import NonlinCorrector, read_extracted_states, petcorr
from simudark_module import simudark_orbvar_function_

nlc = NonlinCorrector()
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
    return (dark*pet + ao)*coadd    

def scia_dark_residuals1(p, data):
    x, y, yerr = data 
    return y - scia_dark_fun1(p, x)

def scia_dark_residuals1e(p, data):
    x, y, yerr = data 
    return (y - scia_dark_fun1(p, x)) / yerr

def extract_05_10_dark_states(orbit):
    # orbit 1
    orbrange = [orbit,orbit]
    states = read_extracted_states(orbrange, 8, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states['mtbl']
    print(states['readoutMean'].shape)
    readouts = states['readoutMean']
    noise = states['readoutNoise']
    n_exec1_8 = readouts.shape[0]
    for idx_exec in range(n_exec1_8):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    pet8 = states['pet'][0]
    coadd8 = states['coadd'][0]
    readouts8_ = readouts[:,7*1024:8*1024] #/ coadd (was already done)
    noise8_ = noise[:,7*1024:8*1024] / numpy.sqrt(coadd8)
    state8_phases = state_mtbl['orbitPhase'][:]

    states63 = read_extracted_states(orbrange, 63, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states63['mtbl']
    print(states['readoutMean'].shape)
    readouts = states63['readoutMean']
    noise = states63['readoutNoise']
    n_exec1_63 = readouts.shape[0]
    for idx_exec in range(n_exec1_63):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    pet63 = states63['pet'][0]
    print('pet63=',pet63)
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
    states = read_extracted_states(orbrange, 8, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states['mtbl']
    print(states['readoutMean'].shape)
    readouts = states['readoutMean']
    noise = states['readoutNoise']
    n_exec2_8 = readouts.shape[0]
    for idx_exec in range(n_exec2_8):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    readouts2_8 = readouts[:,7*1024:8*1024] #/ coadd (was already done)
    noise2_8 = noise[:,7*1024:8*1024] / numpy.sqrt(coadd8)
    state_phases2 = state_mtbl['orbitPhase'][:]

    states63 = read_extracted_states(orbrange, 63, fname, readoutMean=True, readoutNoise=True)
    state_mtbl = states63['mtbl']
    print(states['readoutMean'].shape)
    readouts = states63['readoutMean']
    noise = states63['readoutNoise']
    n_exec2_63 = readouts.shape[0]
    for idx_exec in range(n_exec2_63):
        readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
    pet63 = states63['pet'][0]
    print('pet63=',pet63)
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
    coadd = numpy.concatenate((numpy.zeros(n_exec1_8)+pet8, numpy.zeros(n_exec1_63)+coadd63, numpy.zeros(n_exec2_8)+coadd8, numpy.zeros(n_exec2_63)+coadd63))

    return n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas

# fit dark model to two neighbouring monthly calibration orbits
def fit_monthly(orbit):

    print("ready")

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas = extract_05_10_dark_states(orbit)

    #
    # fit it
    #

    print("steady")

    x = all_state_phases, pet, coadd
    print(all_state_phases.shape)
    print(pet.shape)
    print(pet)
    print(coadd.shape)
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
    amp1_0 = .005
    trend0 = 0
    amp2_0 = 0.1
    phase_offset1 = 0
    phase_offset2 = 0
    p0 = numpy.array([ao0, dc0, amp1_0, trend0, phase_offset1]) #, phase_offset2

    print("go!")

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
        print(pixnr, fitobj.params[4] % 1, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_phases[pixnr] = fitobj.params[4]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()
        #else:
        #   print("Optimal parameters: ", fitobj.params)

    channel_phase = numpy.median(res_phases[numpy.where(statuses >= 0) and numpy.where(res_phases != -3)] % 1)
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

    print("done")

    return channel_phase, aos, lcs, amps, trends

# fit dark model to an orbits with normal (eclipse) dark states
# use parameters computed from nearest calibration orbits
def fit_eclipse_orbit(orbit, aos, lcs, amps, channel_phaseshift):

    print("ready")

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas = extract_05_10_dark_states(orbit)

    #
    # fit it
    #

    print("steady")

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

    print("go!")

    # ..and fit
    n_done = 0
    statuses = numpy.zeros(n_pix)
    res_phases = numpy.zeros(n_pix)
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
        print(pixnr, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_phases[pixnr] = fitobj.params[4]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()
        #else:
        #   print("Optimal parameters: ", fitobj.params)

    print("done")

    return x, lcs, trends, all_readouts, all_sigmas

def compute_trend(orbit, aos):

    n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas = extract_05_10_dark_states(orbit)

    #


#- main ------------------------------------------------------------------------

channel_phase, aos, lcs, amps, trends = fit_monthly(27022)

print('channel_phase=', channel_phase)
print('aos=', aos)
print('lc=', lcs)
print('amp (lc fraction)=', amps)
print('trend (lc fraction)=', trends)

#plt.cla()
#plt.scatter(numpy.arange(1024), aos, c='b')
#plt.scatter(numpy.arange(1024), lcs, c='g')
#plt.show()

x, lcs, trends, readouts, sigmas = fit_eclipse_orbit(27050, aos, lcs, amps, channel_phase)

pts_per_orbit = 50
n_orbits = 2
total_pts = n_orbits*pts_per_orbit
orbphase = numpy.arange(total_pts)/float(pts_per_orbit)
pet05 = numpy.ones(total_pts)*.5 - petcorr
pet10 = numpy.ones(total_pts)*1.0 - petcorr
coadd = numpy.ones(total_pts)

readout_phases, readout_pets, readout_coadd = x

pixnr = 590
x05 = orbphase, pet05, coadd
x10 = orbphase, pet10, coadd
p = aos[pixnr], lcs[pixnr], amps[pixnr], trends[pixnr], channel_phase
plt.cla()
plt.plot(orbphase, scia_dark_fun1(p, x05))
plt.plot(orbphase, scia_dark_fun1(p, x10))
plt.errorbar(readout_phases, readouts[:,pixnr], yerr=sigmas[:,pixnr], ls='none', marker='o')
plt.show()

