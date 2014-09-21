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

from sciamachy_module import NonlinCorrector, read_extracted_states
from simudark_module import simudark_orbvar_function_


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
    dark += dc * amp1 * cos(2*pi*(orbit_phase + phase_shift1))
    #print(amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2)))
    #dark += amp1 * amp2 * cos(4*pi*(orbit_phase + phase_shift2))
    dark += dc * trend * orbit_phase
    return (dark*pet + ao)*coadd    

def scia_dark_residuals1(p, data):
    x, y, yerr = data 
    return y - scia_dark_fun1(p, x)

def scia_dark_residuals1e(p, data):
    x, y, yerr = data 
    return (y - scia_dark_fun1(p, x)) / yerr

def scia_dark_fun(p, x):
    """
     Procedure scia_dark_fun

     Purpose:
       parametric description of dark signal given a number of input parameters:
       * analog offset
       * pixel exposure time
       * co-adding factor
       * dark current
       * orbital phase
       * parameters describing orbital variation
       * parameters describing orbit-to-orbit trend
       Intended to be used with CURVEFIT routine or equivalent

     Usage:
     f = scia_dark_fun(p, x)
     
     Inputs:
     
     p:    1d numpy array containing the function parameters
           * p[0:1024]: numpy array, analog offset in BU for 1024 pixels
           * p[1024:2048]: numpy array, dark current in BU/s for 1024 pixels
           * p[2048:3072]: numpy array, orbital variation amplitude in BU/s
           * p[3072:4096]: numpy array, linear orbit-to-orbit trend in BU/s
           * p[4096]: scalar, orbital variation first overtone relative amplitude in BU/s
           * p[4097]: scalar, orbital variation fundamental frequency phase offset in range [0-1]
           * p[4098]: scalar, orbital variation first overtone phase offset in range [0-1]
     x:    tuple of independent variables (numpy arrays of length n_x)
           * x[0]: numpy array, orbit phase in range (can be [0..1] but also beyond for trending)
           * x[1]: numpy array, pixel exposure time in s
           * x[2]: numpy array, co-adding factor
           all three input values repeated for additional spectra
     
     Returns:
     f:    function values:
           numpy.matrix of dimensions [1024, n_x]
     
     Author:
     Pieter van der Meer, SRON, 2014
     Adapted from IDL code by:
     Ralph Snel, SRON, 2008
    """

    # extract parameters from p
    n_pix = 1024
    analog_offset      = numpy.matrix(p[0:n_pix])
    dark_current       = numpy.matrix(p[n_pix:2*n_pix])
    orb_var_amplitude1 = p[2*n_pix:3*n_pix]
    orb_trend          = numpy.matrix(p[3*n_pix:4*n_pix])
    orb_var_amplitude2 = p[4*n_pix]
    phase_offset1      = p[4*n_pix+1]
    phase_offset2      = p[4*n_pix+2]
    orb_var_amp1_mat   = numpy.matrix(orb_var_amplitude1)

    n_pix = analog_offset.size #nr of pixels
    aones = numpy.matrix(numpy.ones(n_pix))
    # transform to matrices 
    analog_offset = numpy.matrix(analog_offset)
    dark_current = numpy.matrix(dark_current)
    orb_var_amp1_mat = numpy.matrix(orb_var_amplitude1)
    orb_trend = numpy.matrix(orb_trend)
    orb_var_amp2_mat = orb_var_amp1_mat*orb_var_amplitude2

    # Extract exposure information and orbit phase from x
    orbit_phase, pixel_exposure_time, co_adding_factor = x
    n_x = orbit_phase.size # nr of datapoints
    xones = numpy.matrix(numpy.ones(n_x))

    #print('orbit_phase=', orbit_phase)
    #print('PET=', pixel_exposure_time)
    #print('co_adding_factor=', co_adding_factor)
    #print('n_pix=', n_pix,' n_x=', n_x)
	
    # Generate an array for the modeled block of detector read-outs:
    darks = numpy.zeros([n_pix, n_x])

    #
    # Calculate the dark current:
    #
  
    #print('dark_current.shape=',dark_current.shape)
    #print('len(orb_var_amplitude1)=',len(orb_var_amplitude1))
    #print('shapes=',darks.shape,dark_current.T.shape,xones.shape)

    darks += dark_current.T * xones

    # Add the orbital variation of the dark current:
    # Fundamental frequency:
    orbital_variation1 = orb_var_amp1_mat.T * cos((orbit_phase+phase_offset1)*2*pi)
    # First overtone, amplitude relative to fundamental frequency amplitude:
    orbital_variation2 = orb_var_amp2_mat.T * cos((orbit_phase+phase_offset2)*4*pi)
    # Total orbital variation:
    orbital_variation  = orbital_variation1 + orbital_variation2
    darks += orbital_variation
    # orbit-to-orbit trend
    darks += orb_trend.T * orbit_phase
    #print('darks (no pet, no ao)=', darks)

    # Multiply with pixel exposure time:
    #print('pet.shape=',pixel_exposure_time.shape)
    #print('aones.shape=',aones.shape)
    pet_mat = aones.T * pixel_exposure_time
    #print('pet_mat.shape=',pet_mat.shape)
    #print('darks.shape=',darks.shape)
    darks *= pet_mat

    # Add the analog offset:
    darks += (analog_offset.T * xones)

    # Multiply with co-adding factor:
    #print('co_adding_factor.shape=',co_adding_factor.shape)
    coad_mat = aones.T * numpy.matrix(co_adding_factor) 
    darks *= coad_mat

    return darks.T

def scia_dark_residuals(p, data):
    """
    input:
     data: list of two: 
     - x: independent variables (list of three equally sized numpy arrays)
     - y: dependent variables 
    """
    print("ITER")
    x, y = data 
    return y - scia_dark_fun(p, x)

#- main ------------------------------------------------------------------------

nlc = NonlinCorrector()

orbit = 27022
fname = '/SCIA/SDMF31/sdmf_extract_calib.h5'
n_pix = 1024

print("ready")

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

pet = numpy.concatenate((numpy.zeros(n_exec1_8)+pet8, numpy.zeros(n_exec1_63)+pet63, numpy.zeros(n_exec2_8)+pet8, numpy.zeros(n_exec2_8)+pet63))
coadd = numpy.concatenate((numpy.zeros(n_exec1_8)+pet8, numpy.zeros(n_exec1_63)+coadd63, numpy.zeros(n_exec2_8)+coadd8, numpy.zeros(n_exec2_63)+coadd63))

#
# fit it
#

print("steady")

if True:

    x = all_state_phases, pet, coadd
    print(all_state_phases.shape)
    print(pet.shape)
    print(pet)
    print(coadd.shape)
    #fitobj = kmpfit.Fitter(residuals=scia_dark_residuals1, data=(x, all_readouts[:,pixnr], all_sigmas[:,pixnr]))
    
    aoinfo = dict(fixed=False, limits=[0,10000])
    lcinfo = dict(fixed=False, limits=[-100000.,+100000.])
    amp1info = dict(fixed=False, limits=[-.5,+.5])
    trendinfo = dict(fixed=False, limits=[-.1,+.1])
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
        print(pixnr, fitobj.params[4] % 1, fitobj.status, fitobj.message)
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
        print(pixnr, fitobj.params[4] % 1, fitobj.status, fitobj.message)
        statuses[pixnr] = fitobj.status
        res_phases[pixnr] = fitobj.params[4]
        if (fitobj.status <= 0):
           print('Error message = ', fitobj.message)
           quit()

    #plt.cla()
    #plt.scatter(all_state_phases, all_readouts[:,pixnr])
    #plt.show()

else:

    # prepare parameters and fitter
    n_pix = 1024
    x = all_state_phases, pet, coadd
    fitobj = kmpfit.Fitter(residuals=scia_dark_residuals, data=(x, all_readouts))

    aoinfo = [dict(fixed=False) for x in range(n_pix)]
    lcinfo = [dict(fixed=False) for x in range(n_pix)]
    amp1info = [dict(fixed=False, limits=[-1000.,+1000.]) for x in range(n_pix)]
    trendinfo = [dict(Fixed=False, limits=[-1000.,+1000.]) for x in range(n_pix)]
    amp2info = dict(Fixed=False, limits=[-1.,+1.])
    phase1info = dict(fixed=False, limits=[-3.,+3.])
    phase2info = dict(fixed=False, limits=[-3.,+3.])
    fitobj.parinfo = aoinfo+lcinfo+amp1info+trendinfo+[amp2info,phase1info,phase2info]

    # prepare initial parameters
    ao0 = numpy.zeros(n_pix) + 4000
    dc0 = numpy.zeros(n_pix) + 4000
    amp1_0 = numpy.zeros(n_pix) + 10
    trend0 = numpy.zeros(n_pix)
    amp2_0 = 0.1
    phase_offset1 = 0
    phase_offset2 = 0
    p0 = numpy.concatenate((ao0, dc0, amp1_0, trend0, (amp2_0, phase_offset1, phase_offset2)))

    print("go!")

    # ..and fit
    fitobj.fit(params0=p0)
    if (fitobj.status <= 0):
       print('Error message = ', fitobj.message)
    else:
       print("Optimal parameters: ", fitobj.params)

print("done")

