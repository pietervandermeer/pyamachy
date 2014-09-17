import numpy
import h5py
from scipy.optimize import curve_fit

# read simudark data for specified orbit (only ch8 data!)
def read_simudark(orbit, ao=None, lc=None, amp1=None, sig_ao=None, sig_lc=None, sig_amp1=None):
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

    if sig_ao:
        ds_did = gid['sig_ao']
        dict['sig_ao'] = ds_did[:,metaindx]

    if sig_lc:
        ds_did = gid['sig_lc']
        dict['sig_lc'] = ds_did[:,metaindx]

    if sig_amp1:
        ds_did = gid['sig_amp1']
        dict['sig_amp1'] = ds_did[:,metaindx]

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

# version for scipy's curve_fit()
def simudark_orbvar_function_(phases, pet, amp1, amp2, phase1, phase2):

    exec_count = phases.size
    n_pixels   = amp1.size
    func       = numpy.empty((n_pixels, exec_count))

    pet_amp   = numpy.matrix(pet * amp1)
    amp2      = amp2 * pet_amp
    # TODO: expand amp2 to channel width if not already so?!
    #if n_elements(amp2) eq n_pixels/1024 then $
    #    amp2 = rebin(amp2,n_pixels,/sample)

    # correct fundamental frequency
    phases1 = numpy.matrix(phases+phase1) # orbital phases shifted by fundamental phase shift
    func = pet_amp.T * numpy.cos(2*numpy.pi*phases1)

    # correct 1st harmonic
    phases2 = numpy.matrix(phases+phase2)
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

def fit_simudark(x, y):
    popt, pcov = curve_fit(simudark_orbvar_function_, x, y)
    return
