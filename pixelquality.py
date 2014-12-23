#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module containing pixel quality class.
"""

from __future__ import print_function, division

import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import datetime
import warnings
import scipy.signal
warnings.simplefilter("error") # warnings to errors

from read_statedark_module import sdmf_read_statedark
from darklimb_io_module import sdmf_read_rts_darklimb_record
from vardark import load_varkdark_orbit, read_ch8_darks
from sciamachy import get_darkstateid, NonlinCorrector, read_extracted_states, get_closest_state_exec, petcorr, NoiseModel
import config32

#-- functions ------------------------------------------------------------------

def calculate_light_figure(spec, verbose=False, give_smooth=False):
    """
    Calculate sun or wls figure. 
    Smooth (dark and nonlin corrected) spectrum over the pixelrange and check how much the actual pixel readout deviates.
    Deviation is handled absolute-log: exp(-abs(log(smooth/actual)))
    making 0 worst and 1 best case and 1/figure the factor it deviates. 

    Parameters
    ----------

    spec : numpy array, 1D, float
        spectrum (corrected for dark and non-linearity)
    verbose : bool, optional
        if set, this function prints out the worst pixels sorted by descending severity

    Returns
    -------

    reldev : numpy array, dimensions of `spec', float
        quality figures for every pixel in `spec'
    smooth : numpy array, dimensions of `spec', float, optional
        median-smoothed version of `spec'

    """
    smooth = scipy.signal.medfilt(spec) # smooth with median filter
    smooth[smooth==0] = np.nan # avoid division by 0
    raw_figure = np.abs(spec / smooth) # compute ratio between smooth and raw spectrum
    idx = (np.isfinite(raw_figure)) & (np.abs(raw_figure) != 0)
    reldev = np.empty(spec.shape) + np.nan
    reldev[idx] = np.exp(-np.abs(np.log(raw_figure[idx])))

    if verbose:
        sorted_ = np.sort(reldev)
        sorted_pixnr = np.argsort(reldev)
        for i in range(10): # the 10 worst
            print(sorted_pixnr[i], sorted_[i])

    if give_smooth:
        return reldev, smooth
    else:
        return reldev

def create_figure_dset(f, name, dims=None, data=None, chan_sz=1024):
    """ creates a new figure (float per pixel) array in the hdf5 database """
    print("creating "+name)
    if dims is None and data is None:
        dims = (0,chan_sz)
    dtype = np.float64
    ds = f.create_dataset(name, dims, dtype=dtype, #chunks=(16,self.num_chanpixels), 
                          compression='gzip', compression_opts=3, 
                          maxshape=(None,chan_sz), data=data)
    return ds

def create_mask_dset(f, name, dims=None, data=None, chan_sz=1024):
    """ creates a new mask (bool per pixel) array in the hdf5 database """
    print("creating "+name)
    if dims is None and data is None:
        dims = (0,chan_sz)
    dtype = np.bool
    ds = f.create_dataset(name, dims, dtype=dtype, #chunks=(16,self.num_chanpixels), 
                          compression='gzip', compression_opts=3, 
                          maxshape=(None,chan_sz), data=data)
    return ds

def plot_quality_number(spec, thresh):
    """ plot relative deviation and threshold of given channel of pixels """
    thresh_rel = np.array([thresh, thresh])
    thres_x = [0, spec.size]
    thres_y = 1./thresh_rel
    plt.cla()
    plt.plot(np.arange(spec.size), spec, 'bo', thres_x, thres_y, 'k-')
    plt.show()
    return 

class PixelQuality:

    def __init__(self, sdmf30_compat=False):
        
        #
        # define some useful variables
        #

        self.numchannels = 8
        self.num_chanpixels = 1024
        # boundary between 6 and 6+
        self.boundary = self.num_chanpixels*5+795
        self.num_pixels = self.numchannels*self.num_chanpixels
        self.maxorbits = 100000
        self.sdmf30_compat = sdmf30_compat

        #
        # Parse config file, exit if unsuccessful
        #
        
        fname = 'default3.2.cfg'
        cfg_file = open(fname)
        try:
            self.cfg = config32.load(cfg_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+fname)
            raise

        clussizes1 = [       \
        5,192,355,290,177,5, \
        5,71,778,94,71,5,    \
        10,23,897,89,5,      \
        5,5,909,100,5,       \
        5,5,991,18,5,        \
        10,14,973,17,10,     \
        10,38,940,26,10,     \
        10,1004,10           
        ]
        clusoff = np.cumsum(clussizes1)
        self.clusoff = np.insert(clusoff, 0, 0)

        #
        # allocate memory for flagging thresholds
        #

        self.errorlc_thres    = np.zeros(self.num_chanpixels)
        self.residual_thres   = np.zeros(self.num_chanpixels)
        self.exposuretime_max = np.zeros(self.num_chanpixels)
        self.sig_max          = np.zeros(self.num_chanpixels)

        #
        # expand the thresholds now to arrays with size 1024*8
        # this will speed up the calculation later on
        # because i don't have to use loops
        #

        # treat channel 6+ seperately
        range1k = np.arange(self.num_chanpixels)
        for i in range(self.numchannels+1):

            # channels 1 to 5
            if i <= 4: 
                channelindex = range1k + i * self.num_chanpixels
            # channel 6
            if i == 5: 
                channelindex = np.arange(795) + self.num_chanpixels * 5
            # channel 6+
            if i == 6: 
                channelindex = np.arange(229) + self.boundary
            # channels 7 and 8
            if i >= 7: 
                channelindex = range1k + (i-1)*self.num_chanpixels

        #
        # prepare objects relating to low-level detector calibration and statistics
        #

        self.nlc = NonlinCorrector()
        self.noisemodel = NoiseModel()

        #
        # open noise database
        #

        self.noise_fid = h5py.File("/SCIA/SDMF31/pieter/noise.h5", "r")
        self.ds_noise10 = self.noise_fid["pet1.0/noise"]
        self.ds_orbits10 = self.noise_fid["pet1.0/orbits"]
        self.ds_noise05 = self.noise_fid["pet0.5/noise"]
        self.ds_orbits05 = self.noise_fid["pet0.5/orbits"]
        self.ds_noise0125 = self.noise_fid["pet0.125/noise"]
        self.ds_orbits0125 = self.noise_fid["pet0.125/orbits"]

        return

    def get_noise(self, orbit, pet):
        if pet == 0.125:
            return self.get_noise0125(orbit)
        if pet == 0.5:
            return self.get_noise05(orbit)
        if pet == 1.0:
            return self.get_noise10(orbit)

    def get_noise10(self, orbit):
        """ 
        get calibrated noise for orbit. 
        """
        idx = self.ds_orbits10[:] == orbit
        if np.sum(idx) == 0:
            raise Exception("orbit "+str(orbit)+" is not present in noise database.")
        idx = np.argmax(idx)
        return self.ds_noise10[idx,:]

    def get_noise05(self, orbit):
        """ 
        get calibrated noise for orbit. 
        """
        idx = self.ds_orbits10[:] == orbit
        if np.sum(idx) == 0:
            raise Exception("orbit "+str(orbit)+" is not present in noise database.")
        idx = np.argmax(idx)
        return self.ds_noise10[idx,:]

    def get_noise0125(self, orbit):
        """ 
        get calibrated noise for orbit. 
        """
        idx = self.ds_orbits10[:] == orbit
        if np.sum(idx) == 0:
            raise Exception("orbit "+str(orbit)+" is not present in noise database.")
        idx = np.argmax(idx)
        return self.ds_noise10[idx,:]


    def clip_figure(self, figure):
        nonan_idx = np.isfinite(figure)
        if np.sum(nonan_idx) == 0:
            return figure # nan in nan out
        #print("nonan_idx=",nonan_idx)

        idx = figure[nonan_idx] > 1
        if (np.sum(idx) > 0):
            #print("idx=",idx)
            idx_ = np.where(nonan_idx)[0][idx]  # remap index, numpy doesn't do assignments on double sliced stuff
            #print("idx_=",idx_)
            figure[idx_] = 1

        idx = figure[nonan_idx] < 0
        if (np.sum(idx) > 0):
            #print("idx=",idx)
            idx_ = np.where(nonan_idx)[0][idx]  # remap index, numpy doesn't do assignments on double sliced stuff
            #print("idx_=",idx_)
            figure[idx_] = 0
        return figure

    def load_sdmf30_chi(self, orbit):
        fname = "/SCIA/SDMF30/sdmf_dark.h5"
        fid30 = h5py.File(fname, "r")
        grpname = ""
        ds_orbits = fid30[grpname+"orbitList"]
        orbits = ds_orbits[:]
        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")
        i30 = np.argmax(idx)
        dset_chi = fid30["chiSquareFit"]
        return dset_chi[7*1024:8*1024,i30]

    def calculate(self, orbit, debug=False):
        """ 
        perform the quality computation for specified orbit

        Parameters
        ----------

        orbit : int
            absolute orbit number
        """
        
        #
        # get variables from configuration and user arguments
        #

        w_err = self.cfg['w_err']
        w_res = self.cfg['w_res']
        w_noise = self.cfg['w_noise']
        w_sun = self.cfg['w_sun']
        w_wls = self.cfg['w_wls']
        max_noise = self.cfg['max_noise']

        self.orbit = orbit
        darkids = self.cfg['deadflaggingstates']
        calib_db = self.cfg['db_dir']+self.cfg['extract_fname'] # extract_calib database name

        #
        # load vardark
        #

        fname = self.cfg['db_dir']+self.cfg['dark_fname']
        wave_phases, wave_orbit, darkcurrent, analogoffset, uncertainty = load_varkdark_orbit(orbit, False, give_uncertainty=True, fname=fname)
        # slice out a single phase, this is good enough for our purposes
        darkcurrent = darkcurrent[0,0,:].flatten()
        analogoffset = analogoffset.flatten()

        #
        # load sdmf 3.1 dark fit
        #

        fdark = h5py.File("/SCIA/SDMF31/sdmf_dark.h5", 'r')
        gdark = fdark['DarkFit']
        darkcurrent_dset  = gdark["darkCurrent"]
        darkcurrenterror_dset  = gdark["darkCurrentError"]
        analogoffset_dset = gdark["analogOffset"]
        darkcurrent31      = darkcurrent_dset[orbit-1,7*1024:]
#        darkcurrenterror = darkcurrenterror_dset[orbit-1,:]
        analogoffset31     = analogoffset_dset[orbit-1,7*1024:]
        idx = (np.nan_to_num(analogoffset31) > 60000)
        analogoffset31[idx] = np.nan

        #
        # load dark state data (noise/readout for residual)
        #

        fextract = h5py.File(calib_db, 'r')
        readoutnoise = {}
        readoutmean = {}
        readoutpet = {}
        darkpets = (0.125, 0.5, 1.0) # sun, nadir, nadir
        for pet in darkpets:
            darkid = get_darkstateid(pet, orbit)
            jds, readouts, noise, tdet, pets_, coadd = read_ch8_darks([orbit,orbit], darkid)
            pet_string = str(pet)
            readoutmean[pet_string] = readouts
            readoutnoise[pet_string] = noise #self.get_noise(orbit, pet)
            #print(readoutmean[pet_string].shape)

        #
        # check invalid data points in the dark current calculation
        # if the values are 0, the pixel can not be corrected and is 
        # therefore unusable
        #

        self.invalid_mask = (analogoffset == 0) | (np.isnan(analogoffset)) | (darkcurrent == 0) | (np.isnan(darkcurrent))

        #
        # Check for violation of the absolute dark current values
        # These are expressed as a positive and negative saturation time
        # For negative dark currents the threshold defines that the pixel is 
        # not allowed to saturate the adc at adc=0
        # so from the analog offset towards zero within a certain time
        # For positive dark currents to pixel is not supposed to saturate the 
        # adc from analog offset to its maximum value within a certain exposure 
        # time
        #
        # uses: maxadc: max ADC read-out (BU)
        #       sig_max: max signal per channel excluding AO and LC. over the 
        #                sahara, BU/s
        #
        # the quality number variant uses scaling of low and high end of the range
        # to get a linearly varying number between [0..1]
        #

        self.darkcursat_figure = np.empty((1024,), dtype=np.float32)
        exposuretime_max = 1.0-petcorr # for channel 8
        maxadc = 65535
        sig_max = 6000 # max BU/s over sahara. of course this may depend on sensitivity
        lowbg = -analogoffset / exposuretime_max # low dark signal (can be negative)
        highbg_s = (maxadc-analogoffset-sig_max)/exposuretime_max # highest possible dark signal with high signal
        highbg = (maxadc-analogoffset)/exposuretime_max # highest possible dark signal
        darkcurrent_ = np.nan_to_num(darkcurrent)
        darkcursat = (darkcurrent_ < lowbg) | (darkcurrent_ > highbg) # flags (only used for reference..)
        for pixnr in range(1024):
            if self.invalid_mask[pixnr]:
                # invalid pixels are explicitly put to nan here, because this figure cannot be computed for this pixel
                self.darkcursat_figure[pixnr] = np.nan
                continue
            if darkcurrent[pixnr] > highbg_s[pixnr]:
                q = 1.0 - ( (darkcurrent[pixnr] - highbg_s[pixnr]) / (maxadc - highbg[pixnr]) )
            elif darkcurrent[pixnr] < lowbg[pixnr]:
                q = 1.0 - ( (lowbg[pixnr] - darkcurrent[pixnr]) / lowbg[pixnr] )
            else:
                q = 1.0
            if q < 0:
                q = 0
            elif q > 1.0:
                # in theory, this should not be possible. but better safe than sorry. 
                q = 0
            self.darkcursat_figure[pixnr] = q
            #print(pixnr, self.darkcursat_figure[pixnr], analogoffset[pixnr], darkcurrent[pixnr])

        # let's just take a single exposure time to begin with

        self.noise_figure = self.get_noise10(orbit)

        # compute maximum allowable noise based on physical model of on-ground pixels * scalar
        # TODO: replace cfg variable max_noise with noise_scalar
        pixarr = np.arange(7*1024,8*1024, dtype=np.int)
        petarr = np.zeros(1024) + 1.0 - petcorr
        adc = np.abs(np.nan_to_num(darkcurrent_)) * (1.0-petcorr)
        max_noise = 3.0 * self.noisemodel.compute(pixarr, petarr, adc)

        # scale factor for to approximate effect of sigma clipping in sdmf3.0
        if self.sdmf30_compat:
            max_noise *= 1.6

        self.noise_figure = np.nan_to_num(self.noise_figure.flatten())
        idx = self.noise_figure > max_noise
        self.noise_figure[idx] = max_noise[idx]
        # replace bogus 0 value with nan
        idx = np.abs(self.noise_figure) < 0.0001
        self.noise_figure[idx] = np.nan
        self.noise_figure = (max_noise - self.noise_figure) / max_noise
        if debug:
            for pixnr in range(1024):
                print(pixnr, self.noise_figure[pixnr], "<", max_noise[pixnr])

        #
        # check whether the darkcurrent residuals exceed their error
        #

        numdark = len(darkids)
        tmp     = np.zeros(self.num_chanpixels, dtype=np.float64)
        tmp_count = np.zeros(self.num_chanpixels, dtype=np.float64)
        for pet in darkpets:
            pet_string = str(pet)
            state_string = "State_"+format(darkid, "02d")

            #
            # correct dark measurements with dark fit (sounds way too funny)
            #
            
            pets = np.zeros(1024) + pet
            if readoutmean[pet_string].ndim != 2:
                print("ERROR, wrong nr of dimensions!", readoutmean[pet_string], readoutmean[pet_string].shape)
                self.combined = np.ones(self.num_chanpixels)
                return
            corrmean = readoutmean[pet_string] 
            corrmean -= darkcurrent * pets + analogoffset
            corrnoise = readoutnoise[pet_string]
            #print(corrmean)
            #print(corrmean.shape, corrnoise.shape, pets.shape, darkcurrent.shape, analogoffset.shape)

            i_row = 0
            for meanrow in corrmean:
                noiserow = corrnoise[i_row,:]
                #phase = correcteddata[index[j]].phase
                #if phase > 0.0 and phase < 0.3:
                noiserow = np.nan_to_num(noiserow)
                goodnoise = noiserow > 0
                # if there are invalid noise figures, generate an error!
                if np.sum(goodnoise) == 0:
                    logging.warning('invalid noise(=0) in residual criterion!')
                    invalid_mask[:] = True
                tmp[goodnoise] += abs(meanrow[goodnoise] / noiserow[goodnoise])
                tmp_count += goodnoise
                i_row += 1
                
        self.residual_figure = np.sqrt(tmp.flatten() / 20) # TODO: divisor in config file
        self.residual_figure = self.clip_figure(self.residual_figure)

        #
        # compute dark current error to darkcurrent ratio
        #

        self.dc_err_figure = np.abs(darkcurrent / uncertainty)
        #print(self.dc_err_figure.shape)
        self.dc_err_figure = self.dc_err_figure.flatten() / 5000 # empirical scalar to get range normalized to [0..1], TODO: config file
        #print(self.dc_err_figure.shape)
        self.dc_err_figure = self.clip_figure(self.dc_err_figure)

        #
        # prepare for figures that require light measurements (and short vardark product)
        # TODO: disabled until vardark_short is available. SDMF3.1 dark fit will also work well on the short exposure times
        #

        # fname = self.cfg['db_dir']+"pieter/vardark_short.h5"
        # wave_phases_s, wave_orbit_s, darkcurrent_s, analogoffset_s = load_varkdark_orbit(orbit, False, fname=fname)
        # # slice out a single phase, this is good enough for our purposes
        # darkcurrent_s = darkcurrent_s[0,0,:].flatten()
        # analogoffset_s = analogoffset_s.flatten()

        #
        # compute wls figure
        #

        dictwls = get_closest_state_exec(orbit, 61, calib_db, readoutMean=True)
        if debug:
            print("wls pets = ", dictwls['pet'])
        wls_readout = self.nlc.correct_ch8(dictwls['readoutMean']).flatten()
        wls_readout -= darkcurrent31 * dictwls['pet'] + analogoffset31
        self.wls_reldev, smooth_wls = calculate_light_figure(wls_readout, verbose=debug, give_smooth=True)
        #plot_quality_number(self.wls_reldev, self.cfg["thresh_wls"])

        #
        # compute sun figure
        #

        dictsun = get_closest_state_exec(orbit, 62, calib_db, readoutMean=True)
        if debug:
            print("sun pets = ", dictsun['pet'])
        sun_readout = self.nlc.correct_ch8(dictsun['readoutMean']).flatten()
        sun_readout -= darkcurrent31 * dictsun['pet'] + analogoffset31
        self.sun_reldev, smooth_sun = calculate_light_figure(sun_readout, verbose=debug, give_smooth=True)
        # relative threshold
        #plot_quality_number(self.sun_reldev, self.cfg["thresh_sun"])

        #
        # compute sdmf3.0 chisquare (sdmf3.2 vardark is more robust and hence generates different residuals)
        #

        chisquare = self.load_sdmf30_chi(orbit)
        self.chisquare30_figure = (50 - chisquare) / 50 # taken from SDMF3.0 configuration
        self.chisquare30_figure = self.clip_figure(self.chisquare30_figure)

        #
        # combine the criteria using a weighted sum 
        #

        if self.sdmf30_compat:
            self.combined = self.noise_figure * self.noise_figure
            self.combined = self.combined.flatten()
        else:
            self.combined = np.ones(1024, dtype=np.float)
        self.combined *= np.logical_not(self.invalid_mask)
        self.combined *= self.sun_reldev
        self.combined *= self.wls_reldev
        # if SDMF3.0 compatibility mode is enabled:
        if self.sdmf30_compat:
            # chisquare from SDMF 3.0 dark fit
            self.combined *= self.chisquare30_figure
        else:
            # error and residual from SDMF3.2 dark fit
            self.combined *= w_err*self.dc_err_figure + w_res*self.residual_figure 

        #
        # NaN is bogus => quality = 0, clip to range [0..1], and reverse (in preparation for pixelmask (1:bad, 0:good))
        #

        self.combined_flags = np.nan_to_num(self.combined) < 0.005
        if debug:
            for pixnr in range(1024):
                print(pixnr, self.combined_flags[pixnr], 
                             self.combined[pixnr], 
                             self.invalid_mask[pixnr], 
                             self.noise_figure[pixnr], 
                             self.residual_figure[pixnr], 
                             self.dc_err_figure[pixnr],
                             self.wls_reldev[pixnr], 
                             self.sun_reldev[pixnr])

        return

    def write_ascii(self, directory=None, all_figures=False):
        """ 
        write combined quality figure for this orbit to separate ascii file (<orbit>.txt)
        
        Parameters
        ----------

        directory : string, optional
            name of directory (excluding trailing separator)
        """
        if directory is None:
            directory = self.cfg['db_dir']
        fout = open(directory+"/"+str(self.orbit)+".txt", "w")
        for i in range(self.num_chanpixels):
            if all_figures:
                fout.write(str(self.noise_figure[i])+"\t"+str(self.residual_figure[i])+"\t"+str(self.dc_err_figure[i])+"\t"+str(self.combined[i])+"\n")
            else:
                fout.write(str(self.combined[i])+"\n")
        return

    def write(self, directory=None, fname=None):
        """ 
        write data for this orbit to database 
        
        Parameters
        ----------

        directory : string, optional
            name of directory (excluding trailing separator)
        """
        
        orbit = self.orbit
        if fname is not None:
            db_fname = fname
        else:
            if directory is None:
                directory = self.cfg['db_dir']
            db_fname = directory + "/" + self.cfg['pixelmask_fname']
        
        #
        # if database doesn't exist, then create it.. may take a while, but
        # subsequent writes should be quick
        #

        now = datetime.datetime.now()
        nowstr = now.strftime("%Y-%m-%d %H:%M:%S")

        try:
            open(db_fname)
        except IOError as e:
            print("creating db "+db_fname+"...")
            f = h5py.File(db_fname,'w')
            f.attrs["sdmf30_compat"] = self.sdmf30_compat
            orbits = np.array([orbit], dtype='u2')
            entry_dates = np.array([nowstr]) 
            ds = f.create_dataset("orbits", (0,), dtype='u2', 
                             chunks=(1024,), compression='gzip', 
                             compression_opts=3, maxshape=(None,))
            ds.attrs["long_name"] = np.string_("Absolute orbit number")
            f.create_dataset("entryDate", (0,), dtype='S20', 
                             chunks=(1024,), compression='gzip', 
                             compression_opts=3, maxshape=(None,))
            ds = create_figure_dset(f, "wlsResponse", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("WLS Response quality figure (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that gives relative deviation from expected WLS (White Light Source) response.
                                                 For instance, a 10%% response is represented as 0.1, a 1000%% response is too.
                                                 """)

            ds = create_figure_dset(f, "sunResponse", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Sun Response quality figure (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that gives relative deviation from expected Sun-over-ESM diffuser (state 62) response.
                                                 For instance, a 10%% response is represented as 0.1, a 1000%% response is too.
                                                 """)

            ds = create_figure_dset(f, "saturation", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Saturation quality figure (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if the dark current is fully saturated (1.0), 
                                                 not saturated (0.0), or somewhere inbetween.
                                                 This figure is only nonzero if darks are closer than a few thousand BU removed from saturation.""")

            ds = create_figure_dset(f, "darkError", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Dark correction error quality figure (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if the dark correction error (darks corrected by vardark) 
                                                 is too big (1.0), none (0.0), or somewhere inbetween.""")

            ds = create_figure_dset(f, "darkResidual", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Dark fit residual quality figure (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if the fit residual is too big (1.0), none (0.0), or somewhere inbetween.""")

            ds = create_figure_dset(f, "noise", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Noise criterion from SDMF3.0 (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Quality figure [0.0..1.0] that indicates if a pixel exceeds its expected noise 
                                                 (from beginning of mission) by a large margin. 
                                                 This is SDMF3.0 data, which may be used for backward compatibility.""")

            ds = create_figure_dset(f, "combined", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Combined figure (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("Figure [0.0..1.0] that combines all the other figures to estimate a total quality for the pixel (channel 8).")

            ds = create_figure_dset(f, "chisquare3.0", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("Dark Chi^2 from SDMF3.0 (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Chi-squared from the dark signal fit for channel 8. 
                                                 This is SDMF3.0 data, which may be used for backward compatibility.""")

            ds = create_mask_dset(f, "invalid", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("invalid flag (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("Boolean flag that indicates if the dark fit has succeeded (1) or failed (0).")

            ds = create_mask_dset(f, "combinedFlag", chan_sz=self.num_chanpixels)
            ds.attrs["long_name"] = np.string_("combined flag (channel 8)")
            ds.attrs["units"] = np.string_("-")
            ds.attrs["description"] = np.string_("""Boolean flag that indicates a pixel as good (0), or bad (1). 
                                                 It combines all criteria (= `combined' figure thesholded).""")

            f.close()
            print('created db')
        
        #
        # store record in database
        #

        if h5py.h5f.is_hdf5(db_fname):
            f = h5py.File(db_fname)

            if self.sdmf30_compat != f.attrs["sdmf30_compat"]:
                raise Exception("db file has different SDMF3.0 compatibility than PixelQuality object!")
            
            # check if orbit already present..
            dset = f['orbits']
            idx = dset[:] == orbit
            self.n_write = dset.size
            if np.sum(idx) > 0:
                # replace
                idx = np.where(idx)[0][0]
            else:
                # append: resize and write.. 
                self.n_write += 1
                idx = self.n_write-1

            dset.resize((self.n_write,))
            dset[idx] = orbit
            dset = f['entryDate']
            dset.resize((self.n_write,))
            dset[idx] = nowstr
            dset = f['combined']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.combined
            dset = f['darkError']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.dc_err_figure
            dset = f['darkResidual']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.residual_figure
            dset = f['sunResponse']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.sun_reldev
            dset = f['wlsResponse']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.wls_reldev
            dset = f['noise']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.noise_figure
            dset = f['invalid']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.invalid_mask
            dset = f['combinedFlag']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.combined_flags
            dset = f['saturation']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.darkcursat_figure
            dset = f['chisquare3.0']
            dset.resize((self.n_write, self.num_chanpixels))
            dset[idx,:] = self.chisquare30_figure

            f.close()
        else:
            logging.exception("failed to open database %s" % db_fname)
            raise

        return

    def calculate_rts_rank(self, orbit, channel=6, test=False, 
                           pieter_flagging=True, christian_flagging=False, 
                           timedomain_flagging=False, eclipse_flagging=False):
        """
        generate pixelmasks for specified orbit list and channel
        
        Parameters
        ----------
        orbit : int
            absolute orbit number
        channel : int, optional
            6 or 8, 6 implies channel 6+, default = 6
        test : boolean, optional
            set this to plot bar graph instead of writing to db.
        pieter_flagging : boolean, optional
            set this to use time domain flagging, "pieter" style 
        christian_flagging : boolean, optional
            set this to use time domain flagging, "christian" style 
        timedomain_flagging : boolean, optional
            set this to use both pieter and christian flagging
        eclipse_flagging : boolean, optional
            set this to use histogram eclipse data  
        
        Returns
        -------
        combined RTS mask
        
        
        Notes
        ----- 
        - time domain flagging is very strict
        - eclipse flagging involves 3 orbits (!) because this is required to get
          decent statistics. this is probably a bit too strict. 
        """

        #
        # handle user arguments
        #

        if timedomain_flagging:
            pieter_flagging    = True
            christian_flagging = False

        #
        # load eclipse data
        #

        if channel == 6:
            db_name = self.cfg['db_dir']+self.cfg['statedarkch6p_fname']
        if channel == 8:
            db_name = self.cfg['db_dir']+self.cfg['statedarkch8_fname']
        #print 'db_name=', db_name
        data = sdmf_read_statedark([orbit-1,orbit+1], 46, db_name, npeaks=True)
        status = data['status']

        # TODO: exception? or not?
        if status < 0:
            print('error calling sdmf_read_statedark_!')
            return

        ec_npeaks = data['Npeaks']
        statedark_mtbl = data['mtbl']

        # asuming absorbit is the first element..
        n_eff_orbits = statedark_mtbl.size
        # inefficient, but i don't know a way to slice it efficiently..
        ec_orbits = np.zeros(n_eff_orbits, dtype=int)
        for i in range(n_eff_orbits):
            ec_orbits[i] = statedark_mtbl[i][0]
        #print 'ec_orbits=', ec_orbits

        if test:
            print('got eclipse rts data')

        #
        # load dark limb data
        #

        if channel == 6:
            calib_db = self.cfg['db_dir'] + self.cfg['darklimbch6_fname']
        if channel == 8:
            calib_db = self.cfg['db_dir'] + self.cfg['darklimbch8_fname']
        data = sdmf_read_rts_darklimb_record(orbit, calib_db, jumpratios=True, 
            peakratios=True, rtslevel=True, execmean=True, execstddev=True)
            
        #print data
        status = data['status']
        if status < 0:
            print('error calling sdmf_read_rts_darklimb_record ['+str(status)+']!')
            return

        dl_rtslevel = data['rtslevel']
        dl_stds = data['execstddev']
        dl_means = data['execmean']
        peakratios = data['peakratios']
        jumpratios = data['jumpratios']

        if test:
            print('got darklimb rts data')

        #
        # process 2D arrays (orbits x pixels)
        #

        time_method_pieter    = np.zeros(1024, dtype=np.byte)
        time_method_christian = np.zeros(1024, dtype=np.byte)

        #
        # compute "RTS ranking" of darklimb histogram, based on empirical criteria
        #
        #print dl_rtslevel.shape
        dl_rtslice = np.transpose(dl_rtslevel[1,:,0]) # checks if 2nd mode is active (i.e. multi-mode, i.e. rts)
        eclipse_histo_method = np.zeros(1024, dtype=np.byte)
        sum_ec_peaks =np.sum(ec_npeaks,axis=1)
        eclipse_histo_method[795:1024] = sum_ec_peaks[0:229] > 3
        dl_rts = dl_rtslice > 0
        msk = -1
        darklimb_histo_method = np.array(((dl_rts*msk) & np.transpose((1 + (jumpratios >= 1) * (peakratios < 2)))), dtype=np.byte)

        combined_arr = np.zeros(1024)

        # criteria combination bitfield (currently 4 least significant bits)
        # ------------DPEC
        # D: darklimb (always on), P: pieter, E: eclipse, C: christian (not recommended)
        flags = 8+int(pieter_flagging)*4+int(eclipse_flagging)*2+int(christian_flagging)

        for pix in range(1024):

            #
            # find possible rts in time domain (pieter's method): 
            # - too noisy executions
            # - step from one execution to the next  
            # 

            stds  = dl_stds[0,pix,:]
            means = dl_means[0,pix,:]
            #stds  = np.reshape(dl_stds[0,pix,:],?)
            #means = np.reshape(dl_means[0,pix,:],?)

            # remove trailing zeroes
            idx = np.where(means == 0)
            if idx[0].size > 1:
                stds  = stds[0:idx[0]]
                means = means[0:idx[0]]
            else:
                if idx[0] >= 0:
                    continue

            dl_med_std = np.median(stds)
            dl_diffs   = (means-np.roll(means,1))[1:]
            idx_step   = np.where(np.abs(dl_diffs) > 2*dl_med_std)
            idx_noise  = np.where(stds > 2*dl_med_std)
            time_method_pieter[pix] = (idx_step[0].size > 0) or (idx_noise[0].size > 0)

            #
            # find possible rts in time domain (christian's method): 
            # - extreme range 
            # (lower and upper bounds check not possible because not dark corrected!) 
            #

            minmean   = np.min(means)
            maxmean   = np.max(means)
            rangemean = maxmean-minmean
            time_method_christian[pix] = rangemean > 35

        #
        # make slices for this orbit for all the criteria, make a combination and save it all. 
        #

        darklimb_histo_slice = darklimb_histo_method
        eclipse_histo_slice  = eclipse_histo_method
        pieter_slice         = time_method_pieter
        christian_slice      = time_method_christian
        #darklimb_histo_slice = reform(darklimb_histo_method)
        #eclipse_histo_slice  = reform(eclipse_histo_method)
        #pieter_slice         = reform(time_method_pieter)
        #christian_slice      = reform(time_method_christian)
        combined = (darklimb_histo_slice > 0)
        if eclipse_flagging:
            combined = combined | eclipse_histo_slice
        if pieter_flagging:
            combined = combined | pieter_slice
        if christian_flagging:
            combined = combined | christian_slice 
        #print 'pieter_flagging=', pieter_flagging
        #print 'pieter_slice=', pieter_slice
        #print 'darklimb_histo_slice=', darklimb_histo_slice
        #print 'combined=', combined
        #combined_arr[:,i] = combined
        # output: flags, darklimb_histo_slice, eclipse_histo_slice, pieter_slice, christian_slice, combined
        return combined

#- main code -------------------------------------------------------------------

if __name__ == '__main__':
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)
    print("Pixelmask unit test:")
    p = PixelQuality(sdmf30_compat=False)
    print("initialised.")

#    for orbit in range(28631,53200):
#    for orbit in range(4151,53200):
    for orbit in range(42999,43001):
        p.calculate(orbit)
        # try:
        #     p.calculate(orbit)
        # except:
        #    logging.warning("calculation failed for orbit %d!" % orbit)
        #    continue

        print("orbit %d computed" % orbit)
        p.write(directory=".")
        print("Pixel mask data for orbit %d written to db." % orbit)

    #a = p.clip_figure(np.array([np.nan, 0, 0.5, 1.0, 1.5]))
    #print(a)
