#!/usr/bin/env python

# (c) SRON - Netherlands Institute for Space Research (2014).
# All Rights Reserved.
# This software is distributed under the BSD 2-clause license.

#
# Code layout:
#  SECTION VERSION
#  - define various versions of the S/W
#  SECTION AUXILIARY CKD
#  - define functions to read and pre-process CKD required for the calibration
#  - defined functions: get_h5_RadSensMoni
#  SECTION READ DATA
#  - define class SDMFreadSun to read the Sun measurements
#  SECTION CALIBRATE DATA
#  - define class 'SMRcalib' to calibrate the Sun measurements
#  SECTION WRITE DATA
#  - define class 'SMRdb' to write calibrated spectra and meta-data to database
#  SECTION DISPLAY DATA
#  - define class 'SMRshow' to display Sun measurements
#  SECTION ARGPARSE
#  - define function 'handleCmdParams' to obtain command-line parameters
#  SECTION MAIN
#  - code to be run as a standalone program
#
from __future__ import print_function
from __future__ import division

import sys
import argparse
import os.path
import warnings
import numpy as np
import numpy.ma as ma
import h5py
import ConfigParser

from sciamachy_module import petcorr
from scia_dark_functions import scia_dark_fun1
from envisat import PhaseConverter
import config32

#-------------------------SECTION VERSION-----------------------------------
_swVersion = {'major': 0,
              'minor': 8,
              'revision' : 4}
_calibVersion = {'major': 1,
                 'minor': 1,
                 'revision' : 0}
_dbVersion = {'major': 0,
              'minor': 9,
              'revision' : 3}

#-------------------------SECTION ERROR CLASSES-----------------------------
class dbError(Exception):
    pass

class calibrationError(Exception):
    pass

class readSunInfo(Exception):
    pass

# fixed copy of ma.median() 
def mamedian(a, axis=None, out=None, overwrite_input=False):
    from numpy import sort
    """
    Compute the median along the specified axis.

    Returns the median of the array elements.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (None) is
        to compute the median along a flattened version of the array.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output
        but the type will be cast if necessary.
    overwrite_input : bool, optional
        If True, then allow use of memory of input array (a) for
        calculations. The input array will be modified by the call to
        median. This will save memory when you do not need to preserve
        the contents of the input array. Treat the input as undefined,
        but it will probably be fully or partially sorted. Default is
        False. Note that, if `overwrite_input` is True, and the input
        is not already an `ndarray`, an error will be raised.

    Returns
    -------
    median : ndarray
        A new array holding the result is returned unless out is
        specified, in which case a reference to out is returned.
        Return data-type is `float64` for integers and floats smaller than
        `float64`, or the input data-type, otherwise.

    See Also
    --------
    mean

    Notes
    -----
    Given a vector ``V`` with ``N`` non masked values, the median of ``V``
    is the middle value of a sorted copy of ``V`` (``Vs``) - i.e.
    ``Vs[(N-1)/2]``, when ``N`` is odd, or ``{Vs[N/2 - 1] + Vs[N/2]}/2``
    when ``N`` is even.

    Examples
    --------
    >>> x = np.ma.array(np.arange(8), mask=[0]*4 + [1]*4)
    >>> np.ma.extras.median(x)
    1.5

    >>> x = np.ma.array(np.arange(10).reshape(2, 5), mask=[0]*6 + [1]*4)
    >>> np.ma.extras.median(x)
    2.5
    >>> np.ma.extras.median(x, axis=-1, overwrite_input=True)
    masked_array(data = [ 2.  5.],
                 mask = False,
           fill_value = 1e+20)

    """
    if not hasattr(a, 'mask') or np.count_nonzero(a.mask) == 0:
        return masked_array(np.median(a, axis=axis, out=out,
                                  overwrite_input=overwrite_input), copy=False)
    if overwrite_input:
        if axis is None:
            asorted = a.ravel()
            asorted.sort()
        else:
            a.sort(axis=axis)
            asorted = a
    else:
        asorted = sort(a, axis=axis)
    if axis is None:
        axis = 0
    elif axis < 0:
        axis += a.ndim

    # print("asorted.shape[axis]=", asorted.shape[axis])
    # for i in np.arange(238):
    #     print(asorted[i,0])
    counts = asorted.shape[axis] - (asorted.mask).sum(axis=axis)
    h = counts // 2
    # create indexing mesh grid for all but reduced axis
    axes_grid = [np.arange(x) for i, x in enumerate(asorted.shape) if i != axis]
#    print("axes_grid=", axes_grid)
    ind = np.meshgrid(*axes_grid, sparse=True, indexing='ij')
    # insert indices of low and high median
    ind.insert(axis, h - 1)
#    print("ind=", len(ind), ind[0].shape, ind[1].shape, ind)
    low = asorted[ind]
#    print("low=", low)
    ind[axis] = h
#    print("ind=", len(ind), ind[0].shape, ind[1].shape, ind)
    high = asorted[ind]
#    print("high=", high)
    # duplicate high if odd number of elements so mean does nothing
    odd = counts % 2 == 1
    if asorted.ndim == 1:
        if odd:
            low = high
    else:
        # print(odd)
        # print(odd.shape)
        # print(low.mask)
        # print(high.mask)
        # print(low.shape,high.shape)
        low.data[odd] = high.data[odd] # was low[odd] = high[odd], caused a crash in python2.7
    return np.ma.mean([low, high], axis=0, out=out)

#-------------------------SECTION AUXILIARY CKD-----------------------------
def get_h5_RadSensMoni( NDF=True, debug=False ):
    fmtRSPM  ='%-dfloat32, %-dfloat32, 8192float32, (%-d,%-d,8192)float64'
    nameRSPM = ('ang_ele', 'ang_asm', 'wvlen', 'sensitivity')

    fid = h5py.File( '/SCIA/share/nadc_tools/key_radsens.h5', 'r' )
    grp = fid['PPG0']
    dset = grp['PPG0']
    ppg0 = dset[:]

    grp = fid['ABS_RAD']
    dset = grp['Axis 1  Wavelength']
    ref_wl = dset[:]        # regrid all data to ref_wl instead of Sun spectrum
    dset = grp['ABS_RAD']
    abs_rad = dset[:].astype('float64')
    abs_rad /= (5.035e8 * ref_wl * ppg0)
    del(ppg0)

    grp = fid['OBM_s_p']
    dset = grp['Axis 1  Wavelength']
    obm_s_p_wl = dset[:]
    dset = grp['OBM_s_p']
    tmp = dset[:].astype('float64')
    obm_s_p = np.empty_like( tmp )
    for nc in range(8):
        i_mn = nc * 1024
        i_mx = i_mn + 1024
        obm_s_p[i_mn:i_mx] = np.interp( ref_wl[i_mn:i_mx], 
                                        obm_s_p_wl[i_mn:i_mx], 
                                        tmp[i_mn:i_mx] )
    del(tmp)

    grp = fid['ELEV_p']
    dset = grp['Axis 1  wavelength']
    elev_p_wl = dset[:]
    dset = grp['Axis 2  elevation angle']
    elev_p_angle = dset[:]
    dset = grp['ELEV_p']
    elev_p = dset[2,:]

    grp = fid['ELEV_s']
    dset = grp['Axis 1  wavelength']
    elev_s_wl = dset[:]
    dset = grp['Axis 2  elevation angle']
    elev_s_angle = dset[:]
    dset = grp['ELEV_s']
    elev_s = dset[2,:]

    elev_p_a0 = np.interp( ref_wl, elev_p_wl, elev_p.astype('float64') )
    elev_s_a0 = np.interp( ref_wl, elev_s_wl, elev_s.astype('float64') )
    elev_a0 = obm_s_p * elev_s_a0 + elev_p_a0

    if NDF:
        grp = fid['NDF']
        dset = grp['Axis 1  Wavelength']
        ndf_wl = dset[:]
        dset = grp['NDF']
        tmp = dset[:].astype('float64')
        ndf = np.empty_like( tmp )
        for nc in range(8):
            i_mn = nc * 1024
            i_mx = i_mn + 1024
            ndf[i_mn:i_mx] = np.interp( ref_wl[i_mn:i_mx], 
                                        ndf_wl[i_mn:i_mx], 
                                        tmp[i_mn:i_mx] )
        del(tmp)

        grp = fid['NDF_s_p']
        dset = grp['Axis 1  Wavelength']
        ndf_s_p_wl = dset[:]
        dset = grp['NDF_s_p']
        ndf_s_p = dset[:]
        tmp = dset[:].astype('float64')
        ndf_s_p = np.empty_like( tmp )
        for nc in range(8):
            i_mn = nc * 1024
            i_mx = i_mn + 1024
            ndf_s_p[i_mn:i_mx] = np.interp( ref_wl[i_mn:i_mx], 
                                            ndf_s_p_wl[i_mn:i_mx], 
                                            tmp[i_mn:i_mx])
        del(tmp)

        abs_rad *= (2 * ndf / (1 + ndf_s_p))
        obm_s_p *= ndf_s_p

    grp = fid['BRDF_p']
    dset = grp['Axis 1  vacuum wavelength']
    brdf_p_wl = dset[:]
    dset = grp['Axis 2  elevation angle']
    brdf_p_ele = dset[:]
    dset = grp['Axis 3  ASM angle']
    brdf_p_asm = dset[:]
    dset = grp['BRDF_p']
    brdf_p = dset[:].astype('float64')

    grp = fid['BRDF_s']
    dset = grp['Axis 1  vacuum wavelength']
    brdf_s_wl = dset[:]
    dset = grp['Axis 2  elevation angle']
    brdf_s_ele = dset[:]
    dset = grp['Axis 3  ASM angle']
    brdf_s_asm = dset[:]
    dset = grp['BRDF_s']
    brdf_s = dset[:].astype('float64')
    fid.close()

    # re-arrange data: switch ang_asm and ang_ele and store both increasing
    dimWave = brdf_p_wl.shape[0]
    dimASM  = brdf_p_asm.shape[0]
    dimELE  = brdf_p_ele.shape[0]
    dimBRDF = dimASM * dimELE
    indx = []
    for ni in range(dimELE):
        indx += list(ni + np.linspace(dimBRDF-dimELE,0,num=dimASM).astype(int))

    # write data to output structure, interpolated to the ele,asm grid
    rspm = np.empty( 1, dtype=fmtRSPM % (dimELE,dimASM,dimELE,dimASM))
    rspm.dtype.names = nameRSPM
    rspm['ang_ele'] = brdf_p_ele
    rspm['ang_asm'] = brdf_p_asm[::-1]
    rspm['wvlen'] = ref_wl

    brdf_p = brdf_p.reshape( dimBRDF, dimWave )
    brdf_s = brdf_s.reshape( dimBRDF, dimWave )
    for ni in range(dimBRDF):
        brdf_p_ni = np.interp( ref_wl, brdf_p_wl, brdf_p[indx[ni],:] )
        brdf_s_ni = np.interp( ref_wl, brdf_s_wl, brdf_s[indx[ni],:] )
        rspm['sensitivity'][0, ni // dimASM, ni % dimASM,:] = \
            abs_rad * (obm_s_p * brdf_s_ni + brdf_p_ni) / elev_a0

    if debug:
        fid = h5py.File( 'scia_key_rspm.h5', 'w' )
        dset = fid.create_dataset( 'rspm', data=rspm )
        fid.close()
    return rspm

#-------------------------SECTION READ DATA---------------------------------
class SDMFextractSun:
    '''
    Read Sciamachy Sun/State 62 SDMF (v3.1) data
    '''
    def __init__(self, state_id=62, sun_db='/SCIA/SDMF31/sdmf_extract_sun.h5'):
        self.sun_db = sun_db
        self.state_id = state_id
        self.numChannels = 8                 # add some constants for SCIA
        self.channelSize = 1024
        self.numPixels = self.numChannels * self.channelSize
        self.smr = None
        self.wvlen = None

        #
        # Parse config file
        #
        
        fname = 'default3.2.cfg'
        cfg_file = open(fname)
        try:
            self.cfg = config32.load(cfg_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+fname)
            raise

        #
        # read cluster definition
        #

        with h5py.File( sun_db, 'r' ) as fid:
            dset = fid['ClusDef']
            self.clusDef = dset[:]

        return

    def selectOrbits(self, orbitRange):
        with h5py.File( self.sun_db, 'r' ) as fid:
            print("sun_db", self.sun_db)
            grp = fid['State_%02d' % self.state_id]
            dset = grp['orbitList']
            orbitList = dset[:]

        if isinstance(orbitRange,int):
            self.metaIndx = np.argmin(abs(orbitList - orbitRange))
        else:
            self.metaIndx = np.where( (orbitList >= orbitRange[0])
                                      & (orbitList <= orbitRange[1]) )[0]
        self.orbitList = orbitList[self.metaIndx]
        if self.orbitList.size == 0:
            print( '* Info: no orbits selected from sdmf_extract_sun.h5' )
            raise readSunInfo

    def readData(self):
        if self.metaIndx is None or self.metaIndx.size == 0:
            raise readSunInfo

        if isinstance(self.metaIndx, np.ndarray):
            metaIndx = self.metaIndx[0]
            self.metaIndx = self.metaIndx[1:]
        else:
            metaIndx = self.metaIndx
            self.metaIndx = None

        fid = h5py.File( self.sun_db, 'r' )
        grp = fid['State_%02d' % self.state_id]

        dset = grp['metaTable']
        self.mtbl = dset[metaIndx]
        self.absOrbit = self.mtbl['absOrbit']
        self.obmTemp = self.mtbl['obmTemp']
        self.detTemp = self.mtbl['detTemp']

        dset = grp['pointing']
        pointing = dset[metaIndx]
        self.julianDay = np.array([x[0] for x in pointing])
        self.asmAngle = np.array([x[1] for x in pointing])
        self.esmAngle = np.array([x[2] for x in pointing])
        self.sunAzim = np.array([x[3] for x in pointing])
        self.sunElev = np.array([x[4] for x in pointing])

        dset = grp['cluster_01']
        self.numSpectra = dset[0].shape[0]

        self.coaddf = np.ones( (self.numPixels,), dtype=np.uint8 )
        self.pet = np.zeros( (self.numPixels,), dtype=float )
        self.spectra = ma.array( np.zeros(
                (self.numSpectra, self.numPixels), dtype='float64'), 
                                 mask=np.zeros(
                (self.numSpectra, self.numPixels), dtype=int),
                                 hard_mask=True )
        for nc in range( self.clusDef.shape[0] ):
            dset = grp['cluster_%02d' % (nc+1)]
            x = self.clusDef[nc]
            self.coaddf[x[2]:x[2]+x[3]] = dset.attrs['coaddf'][0]
            self.pet[x[2]:x[2]+x[3]] = dset.attrs['PET'][0]
            self.spectra[:,x[2]:x[2]+x[3]] = dset[metaIndx].astype('float64')
        fid.close()

#-------------------------SECTION CALIBRATE DATA----------------------------
class SMRcalib:
    '''
    Listing of implemented calibration IDs:
        Mask dead & blinded pixels
        Co-addition division correction
     1. Memory Effect correction
     2. Non-Linearity correction
     3. Background Signal correction
     4. Stray Light correction
     5. Apply fit parameters
     6. Apply mirror model
     7. Radiance correction
     8. Combine scans to Sun Mean Reference (implied for option "db")
    '''
    def __init__(self, darkVersion="sdmf31"):
        self.funclist = (
            'maskDead', 'coaddDivision', 'memoryEffect', 'nonLinearity', 
            'backGround', 'strayLight', 'fitParam', 'mirrorModel',
            'radiance', 'combineSpectra'
            )
        self.funcdict = dict(
            zip((self.funclist), 
                (self.maskDead, self.coaddDivision, 
                 self.memoryEffect, self.nonLinearity, 
                 self.backGround, self.strayLight, 
                 self.fitParam, self.mirrorModel, 
                 self.radiance, self.combineSpectra))
            )

        self.darkVersion = darkVersion

        #
        # Parse config file
        #
        
        fname = 'default3.2.cfg'
        cfg_file = open(fname)
        try:
            self.cfg = config32.load(cfg_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+fname)
            raise

        return

    def maskDead(self, smr, verbose=False):
        '''
        (*) Identifies dead pixels (based on measurements)
            and blinded pixels at start and end of detector array.
        Parameters
        ----------
        None

        Returns
        -------
        mask where the dead/blinded pixels have boolean value True

        Notes
        -----
        None
        '''
        if verbose:
            print( '(*) Perform masking of the dead/bad pixel' )
        smr.errorType = 'F'
        #
        # mask blinded pixels
        #
        smr.blinded = np.empty( (smr.numPixels,), dtype=np.bool )
        smr.blinded[:] = False

        id_list = np.array( list(range(10)) + list(range(1024-10,1024)) )
        smr.blinded[   0+id_list] = True     # channel 1
        smr.blinded[1024+id_list] = True     # channel 2
        smr.blinded[2048+id_list] = True     # channel 3
        smr.blinded[3072+id_list] = True     # channel 4
        smr.blinded[4096+id_list] = True     # channel 5
        smr.blinded[5120+id_list] = True     # channel 6
        smr.blinded[6144+id_list] = True     # channel 7
        smr.blinded[7168+id_list] = True     # channel 8
        #
        # mask dead pixels
        #
        i_masked = smr.spectra.mask.sum()
        smr.spectra = ma.masked_equal( smr.spectra, 0, copy=False )
        if verbose:
            masked = smr.spectra.mask.sum()
            print( '* Info: masked %6.1f pixels/spectrum with zero signal'
                   % ((masked - i_masked) / float(smr.numSpectra)) )
            i_masked = masked 

        smr.spectra = ma.masked_where( (smr.spectra / smr.coaddf) >= 65535.,
                                       smr.spectra, copy=False )
        if verbose:
            masked = smr.spectra.mask.sum()
            print( '* Info: masked %6.1f pixels/spectrum with saturated signal'
                   % ((masked - i_masked) / float(smr.numSpectra)) )
            i_masked = masked 

    def coaddDivision(self, smr, verbose=False):
        '''
        (*) Co-addition division correction divides the data by the number 
        of measurements added in the on-board co-added.
        Parameters
        ----------
        coaddf : number of measurements co-added (C3)
                 dimension should be one or equal to the spatial dimension

        Returns
        -------
        Signal S_out(i,j) corrected for coadding [ADC counts]

        Notes
        -----
        * None
        '''
        if verbose:
            print( '(*) Perform division by co-adding factor' )
        smr.errorType = 'M'

        smr.spectra /= smr.coaddf

    def memoryEffect(self, smr, verbose=False):
        '''
        (1) Memory Effect correction, this is the effect that the current 
        measurement depends on the previous measurement.
        Parameters
        ----------
        c_mem : memory correction parameters [i,j] (C1)

        Returns
        -------
        Signal S_out(i,j) corrected for memory effect

        Notes
        -----
        * approximation first read-out is obviously wrong, but first readout is
        not used anyway...
        * error estimate not implemented
        '''
        if verbose:
            print( '(1) Perform memory correction (Reticon detectors)' )
        smr.errorType = 'A'
        #
        # read memory correction values
        #
        with h5py.File( '/SCIA/share/nadc_tools/MEMcorr.h5', 'r' ) as fid:
            dset = fid['MemTable']
            memtbl = dset[:]
        #
        # apply memory correction
        #
        id_array = np.arange(smr.channelSize)
        for nch in range(5):
            ipx = id_array + nch * smr.channelSize

            coaddf = smr.coaddf[ipx].max()
            sign = np.rint(smr.spectra[0,ipx]).astype('uint16')
            for nspec in range( smr.numSpectra ):
                corr = memtbl[nch, sign]
                sign = np.rint(smr.spectra[nspec,ipx]).astype('uint16')
                if coaddf > 1:
                    for ni in range(1, coaddf):
                        corr += memtbl[nch, sign]
                    corr /= coaddf
                smr.spectra.data[nspec,ipx] -= corr

    def nonLinearity(self, smr, verbose=False):
        '''
        (2) Non-Linearity correction.
        Parameters
        ----------
        c_nlin : non-linearity correction parameters [i,j] (C1/C2)

        Returns
        -------
        Signal S_out(i,j) corrected for non-linearity effect

        Notes
        -----
        * error estimate not implemented
        '''
        if verbose:
            print( '(2) Perform non-linearity correction (Epitaxx detectors)' )
        smr.errorType = 'A'
        #
        # read non-linearity correction values
        #
        with h5py.File( '/SCIA/share/nadc_tools/NLcorr.h5', 'r' ) as fid:
            dset = fid['CurveIndex']
            curveIndex = dset[:]
            dset = fid['nLinTable']
            nlintbl = dset[:]
        #
        # apply non-linearity correction
        #
        id_array = np.arange(smr.channelSize)
        for nch in range(5, 8):
            pixelList = id_array + nch * smr.channelSize
            curves = curveIndex[nch, id_array]
            for nspec in range(smr.numSpectra):
                sign = np.rint(smr.spectra[nspec,pixelList]).astype('uint16')
                smr.spectra.data[nspec,pixelList] -= nlintbl[curves, sign]

    def backGround(self, smr, verbose=False):
        '''
        (3) Background Signal correction, consists of the dark current (DC) 
        and the thermal background (BG_term).
           BS(i,j) = coaddf * pet * (DC + c_ice * QE * BG_therm)  
        Parameters
        ----------
        coaddf  : number of measurements co-added [i,j] (C3)
        pet     : pixel exposure time [i,j] (C3)
        DC      : dark current (C2/C3)
        c_ice   : transmission coefficient of the ice layer (C2)
        QE      : quantum efficiency of the detector (C1/C2)
        BG_term : thermal background, depends on T_det, T_opt, T_cal and T_grating (C3)

        Returns
        -------
        Signal S_out(i,j) corrected for background signal (i,j)

        Notes
        -----
        * error estimate not implemented
        '''
        from math import cos, pi
        import matplotlib.pyplot as plt

        if verbose:
            print( '(3) Perform subtraction of dark signal' )
        smr.errorType = 'A'

        # make a copy to correct Epitaxx PET without modifying the SMR object
        pet = smr.pet.copy()
        pet[5 * smr.channelSize:] -= petcorr

        #
        # read dark correction values
        #
        with h5py.File( '/SCIA/SDMF31/sdmf_dark.h5', 'r' ) as fid:
            grp = fid['/DarkFit']
            dset = grp['metaTable']
            mtbl = dset[:]
            orbit = mtbl['absOrbit']
            # reject these orbits
            orbit[np.where(mtbl['stateCount'] < 3)] = 999999
            metaIndx = np.argmin(orbit - smr.absOrbit)
            dset = grp['analogOffset']
            ao = dset[metaIndx,:]
            dset = grp['darkCurrent']
            lc = dset[metaIndx,:]

        corr = ao + pet * lc

        corrsimu = corr.copy()
        corrvar = corr.copy()

        # simudark, just for reference
        with h5py.File( "/SCIA/SDMF31/BackUp/sdmf_simudark.h5", 'r' ) as fid:
            grp = fid['/ch8']
            dset = grp['orbitList']
            orbitList = dset[:]
            metaIndx = np.argmin(abs(orbitList - smr.absOrbit))
            dset = grp['metaTable']
            mtbl = dset[metaIndx]
            dset = grp['ao']
            ao = dset[:,metaIndx]
            dset = grp['lc']
            lc = dset[:,metaIndx]
            dset = grp['amp1']
            amp1 = dset[:,metaIndx]

            orbvar = cos( 2 * pi * (mtbl['PHASE1'] + smr.mtbl['orbitPhase']) ) \
                + mtbl['AMP2'] * cos( 4 * pi * (mtbl['PHASE2'] 
                                                + smr.mtbl['orbitPhase']) )
    #        orbsig = cos( 2 * pi * (mtbl['PHASE1'] + smr.mtbl['orbitPhase']) ) \
    #            + mtbl['SIG_AMP2'] * cos( 4 * pi * ( mtbl['PHASE2'] 
    #                                                 + smr.mtbl['orbitPhase']) )

            indx = 7 * smr.channelSize + np.arange(smr.channelSize)
            corrsimu[indx] = ao + pet[indx] * (lc + orbvar * amp1)

        fid_dc = h5py.File( '/SCIA/SDMF31/pieter/vardark_long.h5', 'r' )
        orbits_dc = fid_dc['dim_orbit'][:]
        idx_dc_orb = np.argmin(abs(orbits_dc - smr.absOrbit))
        ds_dc = fid_dc["varDark"] # vardark dataset
        phase_sz = fid_dc["dim_phase"].size

        fid_ao = h5py.File( '/SCIA/SDMF31/pieter/interpolated_monthlies_long.h5', 'r' )
        orbits_ao = fid_ao['orbits'][:]
        idx_ao_orb = np.argmin(abs(orbits_ao - smr.absOrbit))
        ds_ao = fid_ao["aos"] # the analog offsets dataset

        jds = smr.mtbl['julianDay']
        orbit_phase, orbits = phase_conv.get_phase(jds, eclipseMode=True, getOrbits=True)
        real_phase = (orbit_phase + orbits) % 1.
        idx_phase = real_phase * phase_sz

        #orbit_phase = smr.mtbl['orbitPhase']
        n_exec = orbit_phase.size
        # get/compute dark current in BU/s, TODO: interpolated (difference would be minimal)
        dc = ds_dc[idx_dc_orb,idx_phase,:]
        # get/compute analog offset in BU
        ao = ds_ao[idx_ao_orb,:]
        # put these guys in channel 8, where they belong
        corrvar[7*1024:] = ao + dc*pet[7*1024:] 

        #print(smr.spectra.shape, corrvar.shape)
        specsimu = smr.spectra - corrsimu
        specvar = smr.spectra - corrvar           
        specnorm = smr.spectra - corr

        #print(corrsimu[7*1024:8*1024], corrvar[7*1024:8*1024], corr[7*1024:8*1024])
        #print(smr.spectra[smr.numSpectra/2, 7*1024:8*1024])

        plt.cla()
        #print(specsimu.shape, specvar.shape, smr.spectra.shape)
        #print(specsimu[smr.numSpectra/2, 7*1024:8*1024])
        plt.plot(np.arange(1024), specsimu[10, 7*1024:8*1024].flatten(), 'bo', label='simudark 3.1')
        plt.plot(np.arange(1024), specvar[10, 7*1024:8*1024].flatten(), 'ro', label='vardark 3.2')
        plt.plot(np.arange(1024), specnorm[10, 7*1024:8*1024].flatten(), 'go', label='generic 3.1')
        # plt.plot(np.arange(1024), corrvar[7*1024:8*1024], 'ro', label='vardark 3.2')
        # plt.plot(np.arange(1024), corrsimu[7*1024:8*1024], 'bo', label='simudark 3.1')
        # plt.plot(np.arange(1024), corr[7*1024:8*1024], 'go', label='generic 3.1')
        plt.legend(loc='best')
        plt.show()

        #
        # just stick with the generic sdmf 3.1 dark correction for now.. 
        # this saves some bad pixels, although the better pixels are probably 1 or 2 BU off.. 
        # if that matters at all on 10000 BU is to be seen. 
        #

        if self.darkVersion == "sdmf31":
            smr.spectra -= corr 
        elif self.darkVersion == "vardark":
            smr.spectra -= corrvar
        elif self.darkVersion == "simudark":
            smr.spectra -= corrsimu
        else:
            raise calibrationError("unknown dark version "+self.darkVersion)

        #
        # masked invalid pixels
        #
        i_masked = smr.spectra.mask.sum()
        tmp = np.array( [~(np.isfinite(corr)),] * smr.numSpectra )
        smr.spectra = ma.masked_where( tmp, smr.spectra, copy=False )
        del(tmp)
        if verbose:
            masked = smr.spectra.mask.sum()
            print( '* Info: masked %6.1f pixels/spectrum with invalid darks'
                   % ((masked - i_masked) / float(smr.numSpectra)) )
            i_masked = masked

        smr.spectra = ma.masked_less( smr.spectra, 1, copy=False )
        if verbose:
            masked = smr.spectra.mask.sum()
            print( '* Info: masked %6.1f pixels/spectrum with too large darks'
                   % ((masked - i_masked) / float(smr.numSpectra)) )
            i_masked = masked

    def strayLight(self, smr, verbose=False):
        '''
        (4) Stray Light correction
        Parameters
        ----------
        M_stray : stray light correction matrix [i,j,x,y] (C1/C2)
        S_in    : rebin(S_in(i,j) + missing spectrum) [x,y]

        Returns
        -------
        Signal S_out(i,j) corrected for stray light

        Notes
        -----
        * error estimate not implemented
        '''
        if verbose:
            print( '(4) Perform subtraction of spectral stray-light' )
        smr.errorType = 'A'
        #
        # read stray-light correction matrix
        #
        with h5py.File( '/SCIA/share/nadc_tools/Straylight.h5', 'r' ) as fid:
            dset = fid['strayMatrix']
            strayMatrix = dset[:]
            dset = fid['strayGhost']
            strayGhost = dset[:]
            dset = fid['grid_in']
            grid_in = dset[:]
            dset = fid['grid_out']
            grid_out = dset[:]

        # calculate derivative of grid_out
        deriv_out = (np.roll(grid_out, -1) - np.roll(grid_out, 1))/2.
        deriv_out[0] = (4 * grid_out[1] - 3 * grid_out[0] - grid_out[2])/2.
        deriv_out[-1] = (3 * grid_out[-1] - 4 * grid_out[-2] + grid_out[-3])/2.

        # obtain lower and upper indices for regridding, per channel
        low_indx = np.zeros( grid_in.shape )
        high_indx = np.zeros( grid_in.shape )
        input_ch = np.floor(grid_in / smr.channelSize)
        for nc in range(smr.numChannels):
            w = (input_ch == nc)
            grid_ch = grid_in[w]

            # trick to omit round-off errors (fast - only integer manipulation)
            ll = np.empty( grid_ch.shape, dtype=np.uint16 )
            ll[1:] = (grid_ch[:-1] + grid_ch[1:]).astype( 'uint16' )
            ll[(ll % 2) == 1] += 1
            ll //= 2
            ll[0] = nc * smr.channelSize
            ul = np.roll( ll, -1 )
            ul[-1] = (nc+1) * smr.channelSize
            low_indx[w] = ll
            high_indx[w] = ul

        # reduce the spectrum, according to grid_in
        # scale_mask: compensates the reduced spectrum for masked read-outs
        # fillings: compensates for the PET w.r.t. 1/16 sec
        spec_r = np.zeros( (smr.numSpectra, grid_in.shape[0]), dtype='float64' )
        for ni in range(grid_in.shape[0]):
            num = ma.count(smr.spectra[:,low_indx[ni]:high_indx[ni]], axis=1)
            num[num < 1] = 1
            scale_mask = num / float(high_indx[ni] - low_indx[ni])
            fillings = 16 * smr.pet[grid_in[ni]]

            spec_r[:,ni] = smr.spectra[:,low_indx[ni]:high_indx[ni]].sum(axis=1)
            spec_r[:,ni] /= (fillings * scale_mask)
#            print( ni, low_indx[ni], high_indx[ni], spec_r[120,ni], 
#                   scale_mask[120], fillings, deriv_out[ni] )

        # reverse channel 2 (using the numpy 'view' method)
        tmp = spec_r[:,(input_ch == 1)]
        spec_r[:,(input_ch == 1)] = tmp[:,::-1]

        # obtain straylight spectrum
        stray_r = np.dot(spec_r, np.transpose(strayMatrix))

        # correct for sampling distance of the output grid
        stray_r /= deriv_out

        # resample straylight spectrum to SCIA spectrum
        stray = np.zeros( (smr.numSpectra, smr.numPixels), dtype='float64' )
        for ns in range(smr.numSpectra):
            stray[ns,:] = np.interp( np.arange(smr.numPixels, dtype='float64' ),
                                     grid_out, stray_r[ns,:] )

        # blank out blinded pixels
        stray[:, smr.blinded] = 0.

        # scale to original PET of spectra
        stray *= ((16 * smr.pet * smr.coaddf) / smr.coaddf.max())

        # reverse channel 2 (using the numpy 'view' method)
        tmp = stray[:,1024:2048]
        tmp[:,:] = tmp[:,::-1]

        # calculate stray-light contribution of the ghosts
        ghosts = np.zeros( (smr.numSpectra, smr.numPixels), dtype='float64' )
        for ng in range( strayGhost.shape[0] ):
            pp = np.arange( strayGhost[ng,3], strayGhost[ng,5], dtype=int )

            pos = np.polyval( strayGhost[ng,2::-1], pp )
            fact = np.polyval( strayGhost[ng,10:6:-1], pp )
            mask = ((pos >= strayGhost[ng,4]) & (pos <= strayGhost[ng,6]) 
                    & (fact > 0))
    
            pos = pos[mask]
            ghost = fact[mask] * smr.spectra[:,pp[mask]]
            if pos[0] > pos[-1]:
                pos = pos[::-1]
                ghost = ghost[:,::-1]
            pixels = np.arange( int(pos[0]), int(pos[-1]), dtype=int )
            for ns in range(smr.numSpectra):
                ghosts[ns,pixels] += np.interp( pixels, pos, ghost[ns,:] )

#        for ni in range(smr.numPixels):
#            print( ni, smr.spectra[120,ni], stray[120,ni], ghosts[120,ni] )

        # blank out blinded pixels
        ghosts[:, smr.blinded] = 0.

        # subtract straylight from spectra
        smr.spectra -= (stray + ghosts)

    def fitParam(self, smr, verbose=False):
        '''
        (5) Correct Sun (state 62) for the intensity change during the scan
        Parameters
        ----------
        * Read fit parameters from sun_fitpars.h5

        Returns
        -------
        * Sun measurements for intensity change during scan (Diffuser/ESM)

        Notes
        -----
        * error estimate not implemented
        '''
        if verbose:
            print( '(5) Perform correction for diffuser effects' )
        smr.errorType = 'M'
        #
        # read fit parameters for the correction of diffuser effects
        #
        with h5py.File( '/SCIA/SDMF31/Auxiliary/sun_fitpars.h5', 'r' ) as fid:
            dset = fid['fit_parameters']
            fit_param = dset[:]

        jd = smr.julianDay.mean() - fit_param[0,0]
        saa = smr.sunAzim[1:-1].mean() - fit_param[1,0]
        sza = smr.sunElev - fit_param[2,0]

        saa_grid = np.arange(fit_param[3,0]) / (fit_param[3,0] - 1) \
            * (fit_param[5,0]-fit_param[4,0]) + fit_param[4,0] - fit_param[1,0]
        sza_grid = np.arange(fit_param[6,0]) / (fit_param[6,0] - 1) \
            * (fit_param[8,0]-fit_param[7,0]) + fit_param[7,0] - fit_param[2,0]

        for ip in range(smr.numPixels):
            if fit_param[0,ip] == 0: continue

            np0 = 19
            np1 = 19 + fit_param[3, ip]
            saa_lut = fit_param[np0:np1, ip]
            np0 += fit_param[3, ip]
            np1 += fit_param[6, ip]
            sza_lut = fit_param[np0:np1, ip]

            saa_val = np.interp( saa, saa_grid, saa_lut )
            sza_val = np.interp( sza, sza_grid, sza_lut )
            slope = saa_val + fit_param[10, ip] + fit_param[11, ip] * jd \
                + fit_param[12, ip] * jd**2
            ymod = (1 + slope * sza) * (1 + sza_val)
            smr.spectra.data[:,ip] /= ymod

    def mirrorModel(self, smr, verbose=False):
        '''(6) Apply Sun mirror model
        Parameters
        ----------
        Returns
        -------
        Notes
        -----
        Implementation TBD
        '''
        if verbose:
            print( '(6) Perform correction for mirror degradation' )
        smr.errorType = 'M'
        pass

    def radiance(self, smr, verbose=False):
        from pynadc.scia import db,lv1
        '''
        (7) Radiance correction, the light which the telescope recieved in 
        the nadir direction, i.e. all the light comming from the Earth in the 
        FoV.
        Parameters
        ----------
        F_t : radiance sensitivity function, obtained from on-ground, 
              updated by monitor [i,j] (C1/C2)
        
        Returns
        -------
        Radiance I(i,j)

        Notes
        -----
        * error estimate not implemented
        '''
        if verbose:
            print( '(7) Perform radiance correction' )
        smr.errorType = 'M'

        # obtain wavelength grid from level 1b product
        orbits = [int(smr.absOrbit)]
        fileList = db.get_product_by_type( prod_type='1',
                                           proc_best=True, 
                                           orbits=orbits )
        if len(fileList) > 0:
            try:
                l1b = lv1.File( fileList[0] )
                l1b.getSRS()
            except scia_l1b.fmtError as e:
                print( e.msg )
                sys.exit(1)

            wv_interp = True
            smr.wvlen = np.array(l1b.srs['wavelength'][0])
            l1b.__del__()
        else:
            wv_interp = False
            smr.wvlen = smr.rspm['wvlen']    # use rspm grid as approximation

        # reverse channel 2 (using the numpy 'view' method)
        tmp = smr.spectra.data[:,1024:2048]
        tmp[:,:] = tmp[:,::-1]

        # make a copy to correct Epitaxx PET without modifying the SMR object
        pet = smr.pet.copy()
        pet[5 * smr.channelSize:] -= 1.18125e-3

        # linear interpolate for obs_ele, obs_asm
        obs_ele = np.array( smr.sunElev )
        obs_asm = -45 - 0.5 * np.array( smr.asmAngle )
        for no in range(smr.numSpectra):
            na = np.argmin( abs(smr.rspm['ang_asm'][0,:] - obs_asm[no]) )
            if obs_asm[no] < smr.rspm['ang_asm'][0, na] and na > 0: na -= 1
            frac_asm = (obs_asm[no] - smr.rspm['ang_asm'][0, na]) \
                / (smr.rspm['ang_asm'][0, na+1] - smr.rspm['ang_asm'][0, na])

            ne = np.argmin( abs(smr.rspm['ang_ele'][0,:] - obs_ele[no]) )
            if obs_ele[no] < smr.rspm['ang_ele'][0, ne] and ne > 0: ne -= 1
            frac_ele = (obs_ele[no] - smr.rspm['ang_ele'][0, ne]) \
                / (smr.rspm['ang_ele'][0, ne+1] - smr.rspm['ang_ele'][0, ne])

            radsensAzi1 = (1-frac_asm) * smr.rspm['sensitivity'][0,ne,na,:] \
                + frac_asm * smr.rspm['sensitivity'][0,ne,na+1,:]
            radsensAzi2 = (1-frac_asm) * smr.rspm['sensitivity'][0,ne+1,na,:] \
                + frac_asm * smr.rspm['sensitivity'][0,ne+1,na+1,:]
            radsens = (1-frac_ele) * radsensAzi1 + frac_ele * radsensAzi2

            if wv_interp:
                tmp = radsens.copy()
                for nc in range(smr.numChannels):
                    i_mn = nc * smr.channelSize
                    i_mx = i_mn + smr.channelSize
                    radsens[i_mn:i_mx] = np.interp( smr.wvlen[i_mn:i_mx],
                                               smr.rspm['wvlen'][0,i_mn:i_mx], 
                                               tmp[i_mn:i_mx] )
                del(tmp)

            smr.spectra.data[no,:] /= (pet * radsens)

    def combineSpectra(self, smr, verbose=False):
        '''
        (8) Calculate Sun Mean Reference Spectrum
        Parameters
        ----------
        Returns
        -------
        * mean of the 238 scans, skip first and last scan
        * variance of the 238 scans, skip first and last scan
        * bad dead pixel mask

        Notes
        -----
        * error estimate not implemented
        '''
        if verbose:
            print( '(8) Return average of the Sun spectra' )
        smr.errorType = 'A'

        # omit first and last spectra
#        smr.smr      = ma.median( smr.spectra[1:-1,:], axis=0 )
#        smr.smr      = mamedian( smr.spectra[1:-1,:], axis=0 )
        smr.smr      = np.mean( smr.spectra[1:-1,:], axis=0 )
        
        smr.smrVar   = np.zeros( smr.numPixels, dtype='float64' )
        smr.smrSlope = np.zeros( smr.numPixels, dtype='float64' )
        smr.smrError = np.zeros( smr.numPixels, dtype='float64' )
        smr.bdpm     = (ma.count_masked( smr.spectra[1:-1,:], axis=0 ) > 0)

        x = np.arange(238)
        A = np.vstack([x, np.ones(238)]).T
        ipList = np.arange(smr.numPixels)[~smr.bdpm]
        for ip in ipList:
            y = smr.spectra.data[1:-1, ip]
            m, c = np.linalg.lstsq(A, y)[0]
            smr.smrVar[ip] = np.var(y - (m * x + c))
            smr.smrSlope[ip] = m

#-------------------------SECTION WRITE DATA--------------------------------
class SMRdb:
    def __init__( self, args=None, db_name='./sdmf_smr.h5',
                  truncate=False, calibration=None, verbose=False, darkVersion="sdmf31" ):
        if args:
            self.db_name  = args.db_name
            self.calibration = args.calibration.copy()
            self.truncate = args.truncate
            self.verbose  = args.verbose
            self.darkVersion = args.darkVersion
        else:
            self.db_name  = db_name
            self.calibration = calibration
            self.truncate = truncate
            self.verbose  = verbose
            self.darkVersion = darkVersion

    def checkDataBase(self):
        with h5py.File( self.db_name, 'r' ) as fid:
            mystr = ','.join(list(self.calibration.astype('str')))
            if fid.attrs['calibOptions'] != mystr:
                print( 'Fatal:', 'incompatible calibration options' )
                raise dbError('incompatible calibration options')
            if fid.attrs['darkVersion'] != self.darkVersion:
                print( 'Fatal:', 'incompatible dark version', fid.attrs['darkVersion'], "vs", self.darkVersion )
                raise dbError('incompatible dark version')
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

    def fill_mtbl(self, smr):
        from datetime import datetime 

        fmtMTBL  = \
            'float64,a20,uint16,uint16,float32,float32,float32,float32' \
            + ',float32,float32,float32,float32,8float32'
        nameMTBL = ('julianDay','entryDate','absOrbit','quality',
                    'orbitPhase','longitude','latitude','asmAngle',
                    'esmAngle','sunAzim','sunElev','obmTemp',
                    'detTemp')

        self.mtbl = np.empty(1, dtype=fmtMTBL )
        self.mtbl.dtype.names = nameMTBL
        self.mtbl['julianDay'] = smr.mtbl['julianDay']
        self.mtbl['entryDate'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.mtbl['absOrbit'] = smr.absOrbit
        self.mtbl['quality'] = 0
        self.mtbl['orbitPhase'] = smr.mtbl['orbitPhase']
        self.mtbl['longitude'] = smr.mtbl['longitude']
        self.mtbl['latitude'] = smr.mtbl['latitude']
        self.mtbl['asmAngle'] = smr.asmAngle[1:-1].mean()
        self.mtbl['esmAngle'] = smr.esmAngle[1:-1].mean()
        self.mtbl['sunAzim'] = smr.sunAzim[1:-1].mean()
        self.mtbl['sunElev'] = smr.sunElev[1:-1].mean()
        self.mtbl['obmTemp'] = smr.mtbl['obmTemp']
        self.mtbl['detTemp'] = smr.mtbl['detTemp']

    def create(self, smr):
        with h5py.File( self.db_name, 'w', libver='latest' ) as fid:
            ds = fid.create_dataset( 'orbitList', dtype='uint16',
                                     data=smr.absOrbit.reshape(1,), 
                                     maxshape=(None,), chunks=(512,) )
            ds = fid.create_dataset( 'metaTable', 
                                     data=self.mtbl,
                                     chunks=(16384 // self.mtbl.dtype.itemsize,),
                                     shuffle=True, compression='gzip',
                                     compression_opts=1, maxshape=(None,) )
            ds = fid.create_dataset( 'smr', 
                                     data=smr.smr.reshape(1,smr.numPixels), 
                                     maxshape=(None,smr.numPixels), 
                                     chunks=(8,smr.numPixels), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'smrVariance', 
                                     data=smr.smrVar.reshape(1,smr.numPixels), 
                                     maxshape=(None,smr.numPixels), 
                                     chunks=(8,smr.numPixels), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'smrSlope', 
                                     data=smr.smrSlope.reshape(1,smr.numPixels), 
                                     maxshape=(None,smr.numPixels), 
                                     chunks=(8,smr.numPixels), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'smrError', 
                                     data=smr.smrError.reshape(1,smr.numPixels), 
                                     maxshape=(None,smr.numPixels), 
                                     chunks=(8,smr.numPixels), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            ds = fid.create_dataset( 'bdpm', 
                                     data=smr.bdpm.reshape(1,smr.numPixels), 
                                     maxshape=(None,smr.numPixels), 
                                     chunks=(8,smr.numPixels), 
                                     compression='gzip', compression_opts=1,
                                     shuffle=True )
            if smr.wvlen is not None:
                ds = fid.create_dataset( 'wavelength', 
                                         data=smr.wvlen.reshape(1,smr.numPixels), 
                                         maxshape=(None,smr.numPixels), 
                                         chunks=(8,smr.numPixels), 
                                         compression='gzip', compression_opts=1,
                                         shuffle=True )

            # create attributes in the HDF5 root
            mystr = ','.join(list(self.calibration.astype('str')))
            fid.attrs['calibOptions'] = mystr
            fid.attrs['darkVersion'] = self.darkVersion
            fid.attrs['swVersion'] = \
                '%(major)d.%(minor)d.%(revision)d' % _swVersion
            fid.attrs['dbVersion'] = \
                '%(major)d.%(minor)d.%(revision)d' % _dbVersion
            fid.attrs['calibVersion'] = \
                '%(major)d.%(minor)d.%(revision)d' % _calibVersion

    def append(self, smr):
        with h5py.File( self.db_name, 'r+' ) as fid:
            dset = fid['orbitList']         # orbitList
            ax1 = dset.len()
            dset.resize(ax1+1, axis=0)
            dset[ax1] = smr.absOrbit
            orbitList = dset[:]
            dset = fid['metaTable']         # metaTable
            dset.resize(ax1+1, axis=0)
            dset[ax1] = self.mtbl
            dset = fid['smr']               # SMR
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = smr.smr.reshape(1,smr.numPixels)
            dset = fid['smrVariance']       # SMR(variance)
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = smr.smrVar.reshape(1,smr.numPixels)
            dset = fid['smrSlope']       # SMR(slope)
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = smr.smrSlope.reshape(1,smr.numPixels)
            dset = fid['smrError']          # SMR(error)
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = smr.smrError.reshape(1,smr.numPixels)
            dset = fid['bdpm']              # BDPM
            dset.resize(ax1+1, axis=0)
            dset[ax1,:] = smr.bdpm.reshape(1,smr.numPixels)
            if smr.wvlen is not None:
                dset = fid['wavelength']    # waveLength
                dset.resize(ax1+1, axis=0)
                dset[ax1,:] = smr.wvlen.reshape(1,smr.numPixels)

    def rewrite(self, smr):
        with h5py.File( self.db_name, 'r+' ) as fid:
            dset = fid['orbitList']         # orbitList
            orbitList = dset[:]
            ax1 = np.nonzero(orbitList == smr.absOrbit)[0][0]
            dset = fid['metaTable']         # metaTable
            dset[ax1] = self.mtbl
            dset = fid['smr']               # SMR
            dset[ax1,:] = smr.smr.reshape(1,smr.numPixels)
            dset = fid['smrVariance']       # SMR(variance)
            dset[ax1,:] = smr.smrVar.reshape(1,smr.numPixels)
            dset = fid['smrSlope']       # SMR(slope)
            dset[ax1,:] = smr.smrSlope.reshape(1,smr.numPixels)
            dset = fid['smrError']          # SMR(error)
            dset[ax1,:] = smr.smrError.reshape(1,smr.numPixels)
            dset = fid['bdpm']              # BDPM
            dset[ax1,:] = smr.bdpm.reshape(1,smr.numPixels)
            if smr.wvlen is not None:
                dset = fid['wavelength']    # waveLength
                dset[ax1,:] = smr.wvlen.reshape(1,smr.numPixels)

    def store(self, smr, verbose=False):
        self.fill_mtbl( smr )
        if not h5py.is_hdf5( self.db_name ):
            if verbose: print( '* Info: create new database' )
            self.create( smr )
        elif self.truncate: 
            if verbose: print( '* Info: replace database (and start a new)' )
            self.create( smr )
        else:
            self.checkDataBase()
            with h5py.File( self.db_name, 'r' ) as fid:
                dset = fid['orbitList']
                orbitList = dset[:]

            if np.nonzero(orbitList == smr.absOrbit)[0].size == 0:
                if verbose: print( '* Info: append new SMR to database' )
                self.append( smr )
            else:
                if verbose: print( '* Info: overwrite entry in database' )
                self.rewrite( smr )

#-------------------------SECTION DISPLAY DATA------------------------------
class SMRshow:
    def __init__(self, args):
        self.refCalibID = -1
        if args.effect:
            if args.calibration.size == 0:
                print( 'Fatal, need atleast one calibration option' )
                sys.exit(1)
            elif args.calibration.size == 1:
                self.refCalibID = 0
            else:
                self.refCalibID = args.calibration[-2]

    def showPixel(self, smr, pixelID):
        import matplotlib.pyplot as plt

        plt.figure( 1, figsize=(11.69,8.27) )
        plt.title( 'SMR - State 62 - Orbit %d - Pixel %d'
                   % (smr.absOrbit, pixelID) )
        plt.grid( True )
        plt.plot( smr.spectra[:,pixelID], 'b+-' )
        plt.show()
        
    def showSpectrum(self, args, smr):
        import matplotlib.pyplot as plt

        pixelList = []
        for nch in args.channel:
            pixelList += list(range( (nch-1) * smr.channelSize,
                                     nch * smr.channelSize ))

        fig = plt.figure( 1, figsize=(11.69,8.27) )
        fig.text( 0.5, 0.01, 'Temperatures of OBM and Detectors: %6.2f - '
                  % (smr.obmTemp) + np.array_str(smr.detTemp, precision=2),
                  fontsize=9, horizontalalignment='center',
                  verticalalignment='bottom' )
        plt.title( 'SMR - State 62 - Orbit %d' % smr.absOrbit )
        plt.grid( True )

        if smr.wvlen is not None:
            plt.xlim( [smr.wvlen[pixelList[0]], smr.wvlen[pixelList[-1]]] )
            plt.xlabel( 'Wavelength (nm)')

            for nch in args.channel:
                if (nch % 2) == 0:
                    icol = 'b.'
                else:
                    icol = 'r.'
                i_mx = nch * smr.channelSize
                i_mn = i_mx - smr.channelSize
                if smr.smr is None:
                    plt.plot( smr.wvlen[i_mn:i_mx], 
                              smr.spectra[120,i_mn:i_mx], icol )
                else:
                    plt.plot( smr.wvlen[i_mn:i_mx], smr.smr[i_mn:i_mx], icol )
        else:
            plt.xlabel( 'Pixel number')
            plt.xlim( [pixelList[0],pixelList[-1]] )

            if smr.smr is None:
                plt.plot( pixelList, smr.spectra[120,pixelList], 'b.' )
            else:
                plt.plot( pixelList, smr.smr[pixelList], 'b.' )
            plt.ylabel( 'Signal')

        # remove outliers
        tmp = smr.smr[pixelList].compressed()
        tmp = tmp[np.isfinite(tmp)]
        tmp.sort()
        plt.ylim((0,tmp[99*tmp.size/100]))

        plt.show()

    def showEffect(self, args, smr):
        import matplotlib.pyplot as plt

        pixelList = []
        for nch in args.channel:
            pixelList += list(range( (nch-1) * smr.channelSize,
                                     nch * smr.channelSize ))

        fig, axs = plt.subplots( nrows=2, ncols=1, sharex=True, 
                                 figsize=(11.69,8.27) )
        ax = axs[0]
        effect = smr.refSpectra[120,pixelList].copy()
        if smr.errorType == 'M':
            str_effect = '%d / %d' % (args.calibration[-2], args.calibration[-1])
            if smr.smr is None:
                effect /= smr.spectra[120,pixelList]
            else:
                effect /= smr.smr[pixelList]
        elif smr.errorType == 'A':
            str_effect = '%d - %d' % (args.calibration[-2], args.calibration[-1])
            if smr.smr is None:
                effect -= smr.spectra[120,pixelList]
            else:
                effect -= smr.smr[pixelList]
        else:
            pass
        ax.plot( pixelList, effect, 'b.' )
        ax.set_xlim( [pixelList[0],pixelList[-1]] )
        ax.grid( True )
        ax.set_title( 'SMR - calibration: ' + str_effect )

        ax = axs[1]
        if smr.smr is None:
            ax.plot( pixelList, smr.spectra[120,pixelList], 'b.' )
        else:
            ax.plot( pixelList, smr.smr[pixelList], 'b.' )
        ax.set_xlim( [pixelList[0],pixelList[-1]] )
        ax.grid( True )
        ax.set_title( 'SMR - calibration: '
                      + ','.join(str(x) for x in args.calibration )
)
        fig.suptitle( 'SMR - State 62 - Orbit %d'% smr.absOrbit )
        plt.show()
                
    def screen(self, args, smr):
        if args.pixel:
            self.showPixel( smr, args.pixel )
        elif args.effect:
            self.showEffect(args, smr)
        else:
            self.showSpectrum(args, smr)

#-------------------------SECTION ARGPARSE----------------------------------
def handleCmdParams():
    from argparse import ArgumentParser, ArgumentTypeError
    import re

    def parseCalibList( str ):
        '''
        - always perform 'maskDead' and 'coaddDivision'
        - optional are: (1) 'memoryEffect', (2) 'nonLinearity', 
                      (3)'backGround', (4) 'strayLight', (5) 'fitParam', 
                      (6) 'mirrorModel', (7) 'radiance', (8) 'combineSpectra'
        - note: option 'db' implies 'combineSpectra'
        - internally we count from 0 to 9
        '''
        if str.lower() == 'none':
            return np.arange(2, dtype=np.uint8)
        if str.lower() == 'full':
            return np.arange(10, dtype=np.uint8)

        p = re.compile(r'\D')
        if p.match(str) is not None:
            msg = 'Calibration IDs should be digits'
            raise argparse.ArgumentTypeError(msg)
        id_list = np.array(str.split(','), dtype=np.uint8)
        if id_list.min() < 1 or id_list.max() > 8:
            msg = 'Only calibration IDs between 1 and 8 are allowed'
            raise argparse.ArgumentTypeError(msg)

        # add required calibration steps
        return np.concatenate( ([0,1], id_list+1) )
    
    def parseOrbitList( str ):
        msg1 = "'" + str + "' is not a range or number." \
            + " Expected forms like '20000-25000' or '20000'."
        msg2 = "'" + str + "' is not valid orbit number."

        if str.lower() == 'all':
            return None

        m = re.match(r'(\d+)(?:-(\d+))?$', str)
        if not m:
            raise ArgumentTypeError( msg1 )
        v1 = int(m.group(1))
        if m.group(2):
            v2 = int(m.group(2))
            if v1 < 1 or v2 > 100000:
                raise ArgumentTypeError( msg2 )
            return (v1, v2)
        else:
            return v1

    def parseDarkVersion(string):
        if string not in ("sdmf31", "simudark", "vardark"):
            raise ArgumentTypeError("unknown darkVersion "+string)
        return string

    def parseChannelList( str ):
        msg1 = "'" + str + "' is not a range or number." \
            + " Expected forms like '1-5' or '2'."
        msg2 = "'" + str + "' is not valid channel ID."

        m = re.match(r'(\d+)(?:-(\d+))?$', str)
        if not m:
            raise ArgumentTypeError( msg1 )
        v1 = int(m.group(1))
        v2 = int(m.group(2) or v1)
        if v1 < 1 or v2 > 8:
            raise ArgumentTypeError( msg2 )
        return list( range(v1, v2+1) )

    parser = argparse.ArgumentParser( 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description= 'Perform calibration on Sciamachy State 62 data',
        epilog=SMRcalib.__doc__ )

    parser.add_argument( '-v', '--verbose', action='store_true',
                         help='show documentation of object SMRcalib' )
    mystr = '%(major)d.%(minor)d.%(revision)d' % _swVersion
    parser.add_argument( '-V', '--version', action='version', 
                         version='%(prog)s ' + mystr,
                         help='Display the version number then quit')

    # define subparser to display SMR
    subparsers = parser.add_subparsers( title='subcommands',
                                        dest='subparser_name' )
    parser_db = subparsers.add_parser( 'db',
                                    help='options to archive SMR in database' )
    parser_db.add_argument( 'db_name', nargs='?', type=str,
                            default='./sdmf_smr.h5',
                            help='write to hdf5 database' )    
    parser_db.add_argument( '--orbit', type=parseOrbitList,
                            default=None,
                            help='process data from a orbit range' )
    parser_db.add_argument( '-c', '--calibration', default='full', 
                            type=parseCalibList,
            help='calibration IDs (comma separated), or \"none\" or \"full\"' )
    parser_db.add_argument( '--truncate', action='store_true', default=False,
                            help='destroy database and start a new one' )
    parser_db.add_argument( '-P', '--progressbar', action='store_true', 
                            default=False, help='display progress bar' )
    parser_db.add_argument( '-d', '--darkVersion', default='sdmf31', 
                            type=parseDarkVersion,
            help='dark current correction "sdmf31", "vardark", or "simudark"' )
    parser_show = subparsers.add_parser( 'show',
                                         help='options to display SMR' )
    parser_show.add_argument( '--orbit', type=int, default=12005, 
                              help='process Sun spectrum of given orbit' )
    parser_show.add_argument( '-c', '--calibration', default='full', 
                              type=parseCalibList,
            help='calibration IDs (comma separated), or \"none\" or \"full\"' )
    parser_show.add_argument( '--chan', type=parseChannelList,
                              dest='channel', default='1-8',
                              help='show data of one or more channels' )
    parser_show.add_argument( '--pixel', type=int,
                              default=None, help='show data of one pixel' )
    parser_show.add_argument( '--effect', action='store_true', default=False,
                              help='show effect of last calibration step' )
    return parser.parse_args()

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 40 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{}] {:.2%} {}".format( "#"*block + "-"*(barLength-block), progress, status)
    sys.stdout.write(text)
    sys.stdout.flush()

#-------------------------SECTION MAIN--------------------------------------
if __name__ == '__main__':
    phase_conv = PhaseConverter()
    sys.path.append( '/opt/local/EOS/bin' )
    #
    #
    #
    args = handleCmdParams()
    if args.verbose: print( args )
    #
    # show data
    #
    if args.subparser_name == 'show':
        # create object to show figure
        obj_fig = SMRshow(args)

        # create object with calibration routines
        obj_cal = SMRcalib(darkVersion=args.darkVersion)

        # read data
        try:
            smr = SDMFextractSun()
            smr.selectOrbits( args.orbit )
            smr.readData()

            # read required keydata
            if obj_cal.funclist.index('radiance') in args.calibration:
                smr.rspm = get_h5_RadSensMoni()

            # perform calibration
            for calibID in args.calibration:
                obj_cal.funcdict[obj_cal.funclist[calibID]](smr, verbose=args.verbose)
                if calibID == obj_fig.refCalibID:
                    smr.refSpectra = smr.spectra.copy()
        except readSunInfo:
            sys.exit(1)
        except:
            print( "Unexpected error:", sys.exc_info()[0] )
            raise
        else:
            obj_fig.screen(args, smr)
    #
    # store data
    #
    if args.subparser_name == 'db':
        # create object with calibration routines
        obj_cal = SMRcalib(darkVersion=args.darkVersion)

        # make sure we combine the Sun spectra
        indx = obj_cal.funclist.index('combineSpectra')
        if indx not in args.calibration:
            args.calibration = np.append( args.calibration, indx )

        # create database object
        obj_db = SMRdb(args)

        # read data
        smr = SDMFextractSun()
        try:
            smr.selectOrbits( args.orbit )
        except readSunInfo:
            sys.exit(0)

        # read required keydata
        if obj_cal.funclist.index('radiance') in args.calibration:
            smr.rspm = get_h5_RadSensMoni()

        # process all selected orbits
        if args.progressbar:
            p_todo = float(smr.orbitList.size)
            p_done = smr.orbitList.size - smr.metaIndx.size
            update_progress(p_done / p_todo)

        i_orbit = 0
        while True:
            try:
                smr.readData()
                # perform calibration
                for calibID in args.calibration:
                    obj_cal.funcdict[obj_cal.funclist[calibID]](smr, verbose=args.verbose)
                obj_db.store(smr, verbose=args.verbose)

                if args.progressbar:
                    p_done = smr.orbitList.size - smr.metaIndx.size
                    update_progress(p_done / p_todo)

                print("orbit ", smr.orbitList[i_orbit])

            except readSunInfo:
                sys.exit(1)
            except:
                print( "Unexpected error:", sys.exc_info()[0] )
                raise

            i_orbit += 1

