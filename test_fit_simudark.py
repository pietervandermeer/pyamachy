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
from scipy.optimize import curve_fit

from sciamachy_module import NonlinCorrector, read_extracted_states
from simudark_module import simudark_orbvar_function_

#- functions -------------------------------------------------------------------

def scia_dark_fun(x, *a):
    """
     Procedure scia_dark_fun

     Purpose:
       parametric description of dark signal given a
       number of input parameters:
       * analog offset
       * pixel exposure time
       * co-adding factor
       * dark current
       * orbital phase
       * parameters describing orbital variation
       * parameters describing non-linear response (not implemented)
       Intended to be used with CURVEFIT routine or equivalent

     Usage:
     scia_dark_fun, x, a, f
     
     Inputs:
     
     x:    An n-element vector of independent variables
           * x[0]:         orbit phase in range (can be [0..1] but also beyond for trending)
           * x[1]:         pixel exposure time in s
           * x[2]:         co-adding factor
           all three input values repeated for additional spectra
     a:    A vector containing the function parameters
           * a[0:1023]:    analog offset in BU for 1024 pixels
           * a[1024:2047]: dark current in BU/s for 1024 pixels
           * a[2048:3071]: orbital variation amplitude in BU/s
           * a[3072:4095]: linear orbit-to-orbit trend in BU/s
           * a[4096]:      orbital variation first overtone relative amplitude in BU/s
           * a[4097]:      orbital variation fundamental frequency phase offset in range [0-1]
           * a[4098]:      orbital variation first overtone phase offset in range [0-1]
     
     Returns:
     f:    function values:
           * first dimension 1024 pixel values
           * second dimension x.size/3 spectra
     
     Author:
     Pieter van der Meer, SRON, 2014
     Adapted from IDL code by:
     Ralph Snel, SRON, 2008
    """

    # number of pixels
    print('x.size=',x.size)
    print(x)
    print('len(a)=',len(a))
    print(a)
    n_a = ( len(a) - 3 )/4
    aones = numpy.matrix(numpy.ones(n_a))

    # Extract exposure information and orbit phase from x:
    n_x	= x.size/3
    xones = numpy.matrix(numpy.ones(n_x))

    orbit_phase			= x[0:-1:3]
    pixel_exposure_time	= x[1:-1:3]
    co_adding_factor	= x[2:-1:3]

#     print('orbit_phase=', orbit_phase)
#     print('PET=', pixel_exposure_time)
#     print('co_adding_factor=', co_adding_factor)
#     print('n_a=', n_a,' n_x=', n_x)
	
    # Extract the fit parameters from a:
    analog_offset		= numpy.matrix(a[0:n_a-1])
    dark_current		= numpy.matrix(a[n_a:2*n_a-1])
    orb_var_amplitude1	= a[2*n_a:3*n_a-1]
    orb_trend           = numpy.matrix(a[3*n_a:4*n_a-1])
    orb_var_amplitude2	= a[4*n_a]
    phase_offset1       = a[4*n_a+1]
    phase_offset2       = a[4*n_a+2]
    orb_var_amp1_mat    = numpy.matrix(orb_var_amplitude1)
    orb_var_amp2_mat    = numpy.matrix(orb_var_amplitude2)

    # Generate an array for the modeled block of detector read-outs:
    darks = numpy.zeros([n_a, n_x])

    #
    # Calculate the dark current:
    #

    darks += dark_current.T * xones

    # Add the orbital variation of the dark current:
    # Fundamental frequency:
    orbital_variation1 = orb_var_amp1_mat.T * cos((orbit_phase+phase_offset1)*2*pi)
    # First overtone, amplitude relative to fundamental frequency amplitude:
    orbital_variation2 = numpy.matrix(orb_var_amplitude1 * orb_var_amplitude2).T * cos((orbit_phase+phase_offset2)*4*pi)
    # Total orbital variation:
    orbital_variation  = orbital_variation1 + orbital_variation2
    darks += orbital_variation
    # orbit-to-orbit trend
    darks += orb_trend.T * orbit_phase
    #print('darks (no pet, no ao)=', darks)

    # Multiply with pixel exposure time:
    darks *= (pixel_exposure_time * aones).T

    # Add the analog offset:
    darks += (analog_offset * xones).T

    # Multiply with co-adding factor:
    darks *= (co_adding_factor * aones).T

    # Reform the 2-dimensional array to a 1-dimensional vector:
    f = darks.flatten()

    return f

#- main ------------------------------------------------------------------------

nlc = NonlinCorrector()

orbit = 24000
fname = '/SCIA/SDMF31/sdmf_extract_calib.h5'

# orbit 1
orbrange = [orbit,orbit]
states = read_extracted_states(orbrange, 8, fname, readoutMean=True, readoutNoise=True)
state_mtbl = states['mtbl']
print(states['readoutMean'].shape)
readouts = states['readoutMean']
noise = states['readoutNoise']
print(readouts.shape)
for idx_exec in range(readouts.shape[0]):
    readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
pet = states['pet'][0]
coadd = states['coadd'][0]
readouts_ = readouts[:,7*1024:8*1024] #/ coadd (was already done?)
noise_ = noise[:,7*1024:8*1024] / numpy.sqrt(coadd)
state_phases = state_mtbl['orbitPhase'][:]

# orbit 2
orbrange = [orbit+1,orbit+1]
states = read_extracted_states(orbrange, 8, fname, readoutMean=True, readoutNoise=True)
state_mtbl = states['mtbl']
print(states['readoutMean'].shape)
readouts = states['readoutMean']
noise = states['readoutNoise']
print(readouts.shape)
for idx_exec in range(readouts.shape[0]):
    readouts[idx_exec,:] = nlc.correct(readouts[idx_exec,:])
readouts2_ = readouts[:,7*1024:8*1024] #/ coadd (was already done?)
noise2_ = noise[:,7*1024:8*1024] / numpy.sqrt(coadd)
state_phases2 = state_mtbl['orbitPhase'][:]

# combine orbit data
all_state_phases = numpy.concatenate((state_phases, state_phases2))
all_readouts = numpy.concatenate((readouts_, readouts2_)) #.flatten()
all_sigmas = numpy.concatenate((noise_, noise2_)) #.flatten()

x = numpy.empty([3, all_state_phases.size])
x[0,:] = all_state_phases
x[1,:] = pet
x[2,:] = coadd
print(x)
n_pix = 1024
p0 = numpy.zeros(4096+3)
print('x.shape=',x.shape)
print('all_readouts.shape=',all_readouts.shape)
print('all_sigmas.shape=',all_sigmas.shape)
params, pcov = curve_fit(scia_dark_fun, x, all_readouts, sigma=all_sigmas, p0=p0)

