# -*- coding: iso-8859-1 -*-
# module containing pixel quality class.

from __future__ import print_function, division

import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import datetime
import warnings
warnings.simplefilter("error") # warnings to errors

from read_statedark_module import sdmf_read_statedark
from darklimb_io_module import sdmf_read_rts_darklimb_record
from vardark_module import load_varkdark_orbit
from sciamachy_module import get_darkstateid

class PixelQuality:

    def __init__(self):
        
        #
        # define some useful variables
        #

        self.numchannels = 8
        self.num_chanpixels = 1024
        # boundary between 6 and 6+
        self.boundary    = self.num_chanpixels*5+795
        self.num_pixels  = self.numchannels*self.num_chanpixels
        self.maxorbits   = 100000

        #
        # Parse config file, exit if unsuccessful
        #
        
        fname = 'default3.2.cfg'
        cfg_file = open(fname)
        try:
            self.cfg = self.get_config(cfg_file)
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

    def get_config(self, config_file):
        """ 
        load configuration from file 

        Parameters
        ----------

        config_file : string
            name of the configuration file

        Returns
        -------
        get_config : dict
            contains all relevant configuration settings as readily-usable types (int, boolean, float, arrays instead of strings)
        """
        import string

        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['extract_fname'] = parser.get('Global','extract_file')
        dict['dark_fname'] = parser.get('Global','dark_file')
        dict['pixelmask_fname'] = parser.get('Global','pixelmask_file')
        dict['statedarkch6p_fname'] = parser.get('Global','statedarkch6p_file')
        dict['statedarkch8_fname'] = parser.get('Global','statedarkch8_file')
        dict['darklimbch6_fname'] = parser.get('Global','darklimbch6_file')
        dict['darklimbch8_fname'] = parser.get('Global','darklimbch8_file')
        dict['darklimbch8_fname'] = parser.get('Global','darklimbch8_file')
        string = parser.get('Processor','dc_sat_time')
        dict['dc_sat_time'] = [float(s) for s in string.split(',')]
        string = parser.get('Processor','dc_err_thres')
        dict['dc_err_thres'] = [float(s) for s in string.split(',')]
        string = parser.get('Processor','res_thres')
        dict['res_thres'] = [float(s) for s in string.split(',')]
        string = parser.get('Processor','max_signal')
        dict['max_signal'] = [float(s) for s in string.split(',')]
        string = parser.get('Processor','deadflaggingstates')
        dict['deadflaggingstates'] = [int(s) for s in string.split(',')]
        dict['w_err'] = float(parser.get('Processor','w_err'))
        dict['w_res'] = float(parser.get('Processor','w_res'))
        dict['w_noise'] = float(parser.get('Processor','w_noise'))
        return dict

    def clip_figure(self, figure):
        nonan_idx = np.logical_not(np.isnan(figure))
        if np.sum(nonan_idx) == 0:
            return figure # nan in nan out
        #print("nonan_idx=",nonan_idx)

        idx = figure[nonan_idx] > 1
        if (np.sum(idx) > 0):
            #print("idx=",idx)
            idx_ = np.where(nonan_idx)[0][idx]  # remap index, numpy doesn't do assignments on double sliced stuff
            #print("idx_=",idx_)
            figure[idx_] = 1
        return figure

    def calculate(self, orbit):
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
        self.orbit = orbit
        darkids = self.cfg['deadflaggingstates']

        #
        # load vardark
        #

        fname = self.cfg['db_dir']+self.cfg['dark_fname']
        wave_phases, wave_orbit, darkcurrent, analogoffset, uncertainty = load_varkdark_orbit(orbit, False, give_uncertainty=True, fname=fname)
        # slice out a single phase, this is good enough for our purposes
        darkcurrent = darkcurrent[0,0,:].flatten()
        analogoffset = analogoffset.flatten()

        #
        # load dark state data (noise/readout for residual)
        #

        fextract = h5py.File(self.cfg['db_dir']+self.cfg['extract_fname'], 'r')
        readoutnoise = {}
        readoutmean = {}
        readoutpet = {}
        darkpets = (0.125, 0.5, 1.0) # sun, nadir, nadir
        for pet in darkpets:
            darkid = get_darkstateid(pet, orbit)
            state_string = "State_"+format(darkid, "02d")
            pet_string = str(pet)
            gextract = fextract[state_string]
            readoutnoise_dset = gextract["readoutNoise"]
            readoutmean_dset = gextract["readoutMean"]
            # TODO: nonlinearity correction
            # TODO: SDMF3.1 noise is wrong, need to use SDMF3.0 or my own noise product derived from this!
            orbitlist = gextract["orbitList"]
            idx = np.where(orbitlist[:] == orbit)
            if idx[0].size == 0:
                logging.warning("no dark data for orbit %d!" % orbit)
                self.invalid_mask = np.ones(self.num_chanpixels, dtype=np.bool)
                self.combined = np.empty((self.num_chanpixels,), dtype=np.float32)
                self.combined[:] = np.nan
                return
            readoutnoise[pet_string] = readoutnoise_dset[idx[0],:]
            readoutmean[pet_string] = readoutmean_dset[idx[0],:]

        #
        # check invalid data points in the dark current calculation
        # if the values are 0, the pixel can not be corrected and is 
        # therefore unusable
        #

        self.invalid_mask = (analogoffset == 0) | (np.isnan(analogoffset)) | (darkcurrent == 0) | (np.isnan(darkcurrent))

        # let's just take a single exposure time to begin with

        if "1.0" in readoutnoise:
            if readoutnoise["1.0"].ndim != 2:
                print("ERROR, wrong nr of dimension!", readoutnoise["1.0"], readoutnoise["1.0"].shape)
                self.combined = np.ones(self.num_chanpixels)
                return
            else:
                self.noise_figure = np.mean(readoutnoise["1.0"], axis=0)[7*1024:]
        else:
            self.noise_figure = np.ones(self.num_chanpixels) * np.nan
        # replace bogus 0 value with nan
        idx = self.noise_figure == 0
        if np.sum(idx) > 0:
            self.noise_figure[idx] = np.nan
        self.noise_figure = 1024 / self.noise_figure.flatten()

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
                print("ERROR, wrong nr of dimension!", readoutmean[pet_string], readoutmean[pet_string].shape)
                self.combined = np.ones(self.num_chanpixels)
                return
            corrmean = readoutmean[pet_string][:,7*1024:8*1024] 
            corrmean -= darkcurrent * pets + analogoffset
            corrnoise = readoutnoise[pet_string][:,7*1024:8*1024]
            #print(corrmean)
            #print(corrmean.shape, corrnoise.shape, pets.shape, darkcurrent.shape, analogoffset.shape)

            i_row = 0
            for meanrow in corrmean:
                noiserow = corrnoise[i_row,:]
                #phase = correcteddata[index[j]].phase
                #if phase > 0.0 and phase < 0.3:
                #print(noiserow.shape)
                noiserow = np.nan_to_num(noiserow)
                goodnoise = noiserow > 0
                # if there are invalid noise figures, generate an error!
                if np.sum(goodnoise) == 0:
                    logging.warning('invalid noise(=0) in residual criterion!')
                    invalid_mask[:] = True
                #print(goodnoise.shape, meanrow)
                tmp[goodnoise] += abs(meanrow[goodnoise] / noiserow[goodnoise])
                tmp_count += goodnoise
                i_row+=1
                
        self.residual_figure = np.sqrt(tmp.flatten()) / np.sqrt(20)
        self.residual_figure = self.clip_figure(self.residual_figure)

        #
        # compute dark current error to darkcurrent ratio
        #

        self.dc_err_figure = np.abs(darkcurrent / uncertainty)
        print(self.dc_err_figure.shape)
        self.dc_err_figure = self.dc_err_figure.flatten() / 5000 # empirical scalar to get range normalized to [0..1]
        print(self.dc_err_figure.shape)
        self.dc_err_figure = self.clip_figure(self.dc_err_figure)

        #
        # compute wls figure
        #

        #

        #
        # compute sun figure
        #

        #

        #
        # combine the criteria using a weighted sum 
        #

        print(self.dc_err_figure.shape)
        print(self.residual_figure.shape)
        print(self.noise_figure.shape)
        self.combined = w_err*self.dc_err_figure + w_res*self.residual_figure + w_noise*self.noise_figure
        self.combined = self.combined.flatten()
        #print(self.combined, self.combined.shape, self.combined.dtype)
        self.combined *= np.logical_not(self.invalid_mask)
        self.combined = self.combined*w_sun + (1.-w_sun)*self.combined*self.sun_figure
        self.combined = self.combined*w_wls + (1.-w_wls)*self.combined*self.wls_figure

        #
        # NaN is bogus => quality = 0, clip to range [0..1], and reverse (in preparation for pixelmask (1:bad, 0:good))
        #

        self.combined = np.nan_to_num(self.combined)
        idx = self.combined > 1
        if (np.sum(idx) > 0):
            self.combined[idx] = 1
        self.combined = 1 - self.combined

        return        

    def create_figure_dset(self, f, name):
        """ creates a new figure (float per pixel) array in the database """
        print("creating "+name)
        dims = (0,self.num_chanpixels)
        dtype = np.float64
        f.create_dataset(name, dims, dtype=dtype, #chunks=(16,self.num_chanpixels), 
                         compression='gzip', compression_opts=3, 
                         maxshape=(None,self.num_chanpixels))

    def create_mask_dset(self, f, name):
        """ creates a new mask (bool per pixel) array in the database """
        print("creating "+name)
        dims = (0,self.num_chanpixels)
        dtype = np.bool
        f.create_dataset(name, dims, dtype=dtype, #chunks=(16,self.num_chanpixels), 
                         compression='gzip', compression_opts=3, 
                         maxshape=(None,self.num_chanpixels))

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

    def write(self, directory=None):
        """ 
        write data for this orbit to database 
        
        Parameters
        ----------

        directory : string, optional
            name of directory (excluding trailing separator)
        """
        
        orbit = self.orbit
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
            print('creating meta data')
            orbits = np.array([orbit], dtype='u2')
            entry_dates = np.array([nowstr]) 
            f.create_dataset("orbits", (0,), dtype='u2', 
                             chunks=(1024,), compression='gzip', 
                             compression_opts=3, maxshape=(None,))
            f.create_dataset("entryDate", (0,), dtype='S20', 
                             chunks=(1024,), compression='gzip', 
                             compression_opts=3, maxshape=(None,))
            self.create_figure_dset(f, "darkError")
            self.create_figure_dset(f, "darkResidual")
            self.create_figure_dset(f, "noise")
            self.create_figure_dset(f, "combined")
            self.create_mask_dset(f, "invalid")
            f.close()
            print('created db')
        
        #
        # store record in database
        #

        if h5py.h5f.is_hdf5(db_fname):
            f = h5py.File(db_fname)
            
            # check if orbit already present..
            dset = f['orbits']
            idx = dset[:] == orbit
            if np.sum(idx) > 0:
                # replace
                idx = np.where(idx)[0][0]
                dset[idx] = orbit
                dset = f['entryDate']
                dset[idx] = nowstr
                dset = f['combined']
                dset[idx,:] = self.combined
                dset = f['darkError']
                dset[idx,:] = self.dc_err_figure
                dset = f['darkResidual']
                dset[idx,:] = self.residual_figure
                dset = f['noise']
                dset[idx,:] = self.noise_figure
                dset = f['invalid']
                dset[idx,:] = self.invalid_mask
            else:
                # append: resize and write.. orbits dataset first
                self.n_write = dset.size + 1
                dset.resize((self.n_write,))
                dset[self.n_write-1] = orbit

                # now the other data sets
                dset = f['entryDate']
                dset.resize((self.n_write,))
                dset[self.n_write-1] = nowstr
                dset = f['combined']
                dset.resize((self.n_write, self.num_chanpixels))
                dset[self.n_write-1,:] = self.combined
                dset = f['darkError']
                dset.resize((self.n_write, self.num_chanpixels))
                dset[self.n_write-1,:] = self.dc_err_figure
                dset = f['darkResidual']
                dset.resize((self.n_write, self.num_chanpixels))
                dset[self.n_write-1,:] = self.residual_figure
                dset = f['noise']
                dset.resize((self.n_write, self.num_chanpixels))
                dset[self.n_write-1,:] = self.noise_figure
                dset = f['invalid']
                dset.resize((self.n_write, self.num_chanpixels))
                dset[self.n_write-1,:] = self.invalid_mask
            
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
    p = PixelQuality()
    orbit = 10000
    print("initialised.")
    p.calculate(orbit)
    print("orbit %d computed" % orbit)
    p.write(directory=".")
    print("Pixel mask data for orbit %d written to db." % orbit)

    a = p.clip_figure(np.array([np.nan, 0, 0.5, 1.0, 1.5]))
    print(a)