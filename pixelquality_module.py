# -*- coding: iso-8859-1 -*-
# module containing pixel quality class.

from __future__ import print_function, division

import ConfigParser
import numpy as np
import numpy
import matplotlib.pyplot as plt
import h5py
import logging
import datetime
import warnings
warnings.simplefilter("error") # warnings to errors
numpy.set_printoptions(threshold=numpy.nan, precision=4, suppress=True, linewidth=np.nan)

from read_statedark_module import sdmf_read_statedark
from darklimb_io_module import sdmf_read_rts_darklimb_record
from vardark_module import load_varkdark_orbit

class PixelQuality:

    def __init__(self):
        
        #
        # define some useful variables
        #

        self.numchannels = 8
        self.numpixels   = 1024
        # boundary between 6 and 6+
        self.boundary    = self.numpixels*5+795
        self.pixelcount  = self.numchannels*self.numpixels
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
        clusoff = numpy.cumsum(clussizes1)
        self.clusoff = numpy.insert(clusoff, 0, 0)

        #
        # allocate memory for flagging thresholds
        #

        self.errorlc_thres    = numpy.zeros(self.pixelcount)
        self.residual_thres   = numpy.zeros(self.pixelcount)
        self.exposuretime_max = numpy.zeros(self.pixelcount)
        self.sig_max          = numpy.zeros(self.pixelcount)

        #
        # expand the thresholds now to arrays with size 1024*8
        # this will speed up the calculation later on
        # because i don't have to use loops
        #

        # treat channel 6+ seperately
        range1k = numpy.arange(self.numpixels)
        for i in range(self.numchannels+1):

            # channels 1 to 5
            if i <= 4: 
                channelindex = range1k + i * self.numpixels
            # channel 6
            if i == 5: 
                channelindex = numpy.arange(795) + self.numpixels * 5
            # channel 6+
            if i == 6: 
                channelindex = numpy.arange(229) + self.boundary
            # channels 7 and 8
            if i >= 7: 
                channelindex = range1k + (i-1)*self.numpixels

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

        wave_phases, wave_orbit, darkcurrent, analogoffset = load_varkdark_orbit(orbit, False)
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
        for darkid in darkids:
            state_string = "State_"+format(darkid, "02d")
            gextract = fextract[state_string]
            readoutnoise_dset = gextract["readoutNoise"]
            readoutmean_dset = gextract["readoutMean"]
            orbitlist = gextract["orbitList"]
            idx = numpy.where(orbitlist[:] == orbit)
            if idx[0].size == 0:
                logging.warning("no dark data for orbit %d!" % orbit)
                invalid_mask[:] = True
            readoutnoise[state_string] = readoutnoise_dset[idx[0],:]
            readoutmean[state_string] = readoutmean_dset[idx[0],:]

            #
            # get pets for this orbit
            #

            clusconf = gextract["clusConf"]
            # find row based on orbit number
            pets = numpy.zeros(self.pixelcount)
            idx_row = -1
            for row in clusconf:
                idx_row += 1
                if orbit < row[0]:
                    break
            clusconf_ = clusconf[:]
            row = clusconf_[idx_row-1]
            cluspets = row[3]
            for i_clus in range(40):
                i_start = self.clusoff[i_clus]
                i_end = self.clusoff[i_clus+1]
                pets[i_start:i_end] = cluspets[i_clus]
            
            readoutpet[state_string] = pets
            
        self.readoutpet = readoutpet

        #
        # check invalid data points in the dark current calculation
        # if the values are 0, the pixel can not be corrected and is 
        # therefore unusable
        #

        invalid_mask = (analogoffset == 0) | (np.isnan(analogoffset)) | (darkcurrent == 0) | (np.isnan(darkcurrent))

        #
        # check whether the darkcurrent residuals exceed their error
        #

        numdark = len(darkids)
        tmp     = np.zeros(self.numpixels, dtype=np.float64)
        tmp_count = np.zeros(self.numpixels, dtype=np.float64)
        for darkid in darkids:
            state_string = "State_"+format(darkid, "02d")

            #
            # correct dark measurements with dark fit (sounds way too funny)
            #
            
            pets = readoutpet[state_string][7*1024:8*1024]
            corrmean = readoutmean[state_string][:,7*1024:8*1024] 
            corrmean -= darkcurrent * pets + analogoffset
            corrnoise = readoutnoise[state_string][:,7*1024:8*1024]
            #print(corrmean)
            print(corrmean.shape, corrnoise.shape, pets.shape, darkcurrent.shape, analogoffset.shape)

            i_row = 0
            for meanrow in corrmean:
                noiserow = corrnoise[i_row,:]
                #phase = correcteddata[index[j]].phase
                #if phase > 0.0 and phase < 0.3:
                print(noiserow.shape)
                noiserow = np.nan_to_num(noiserow)
                goodnoise = noiserow > 0
                # if there are invalid noise figures, generate an error!
                if np.sum(goodnoise) == 0:
                    logging.warning('invalid noise(=0) in residual criterion!')
                    invalid_mask[:] = True
                print(goodnoise.shape, meanrow)
                tmp[goodnoise] += abs(meanrow[goodnoise] / noiserow[goodnoise])
                tmp_count += goodnoise
                i_row+=1
                
        residual_figure = tmp

        #
        # compute dark current error to darkcurrent ratio
        #

        dc_err_figure = np.abs(darkcurrenterror / darkcurrent)

        #
        # combine the criteria using a weighted sum 
        #

        self.combined = (w_err*dc_err_figure + w_res*residual_figure + w_noise*noise_figure) * (not invalid_mask)

        return        

    def create_figure_dset(self, f, name, orbit):
        """ creates a new mask array in the database """
        print("creating "+name)
        dims = (self.maxorbits,self.pixelcount)
        if (orbit < 1) or (orbit > self.maxorbits):
            raise ValueError('orbit %d out of range [1,%d]' %
                             (orbit,self.maxorbits))
        dtype = numpy.float64
        dat = numpy.zeros(dims, dtype=dtype)
        f.create_dataset(name, dims, dtype=dtype, chunks=(16,self.numpixels), 
                         compression='gzip', compression_opts=3, 
                         maxshape=[None,self.pixelcount], data=dat)

    def write(self):
        """ write data for this orbit to database """
        
        orbit = self.orbit
        meta_dtype = numpy.dtype([('absOrbit', '=u4'), ('entryDate', '=S20')])
        db_fname = self.cfg['db_dir']+self.cfg['pixelmask_fname']
        
        #
        # if database doesn't exist, then create it.. may take a while, but
        # subsequent writes should be quick
        #

        try:
            open(db_fname)
        except IOError as e:
            print("creating db "+db_fname+"...")
            f = h5py.File(db_fname,'w')
            print('creating metaTable')
            metadata = numpy.empty((self.maxorbits), dtype=meta_dtype)
            f.create_dataset("metaTable", (self.maxorbits,), dtype=meta_dtype, 
                             chunks=(1024,), compression='gzip', 
                             compression_opts=3, maxshape=None, data=metadata)
            self.create_figure_dset(f, "darkError", orbit)
            self.create_figure_dset(f, "darkResidual", orbit)
            self.create_figure_dset(f, "noise", orbit)
            self.create_figure_dset(f, "combined", orbit)
            f.close()
            print('created db')

        #
        # modify record in database
        #

        if h5py.h5f.is_hdf5(db_fname):
            f = h5py.File(db_fname)
            
            now = datetime.datetime.now()
            nowstr = now.strftime("%Y-%m-%d %H:%M:%S")
            metadata = numpy.empty((1), dtype=meta_dtype)
            metadata[0] = (orbit, nowstr)
            
            dset = f['metaTable']
            dset[orbit-1] = metadata
            dset = f['RTS']
            dset[orbit-1,:] = self.rts_mask
            dset = f['combined']
            dset[orbit-1,:] = self.combined_mask
            dset = f['darkCurrentError']
            dset[orbit-1,:] = self.dc_err_mask
            dset = f['darkCurrentSat']
            dset[orbit-1,:] = self.darkcursat_mask
            dset = f['residual']
            dset[orbit-1,:] = self.residual_mask
            dset = f['invalid']
            dset[orbit-1,:] = self.invalid_mask
            
            f.close()
        else:
            logging.exception("failed to open database %s" % db_fname)
            raise

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
        pieter_flagging: boolean, optional
            set this to use time domain flagging, "pieter" style 
        christian_flagging: boolean, optional
            set this to use time domain flagging, "christian" style 
        timedomain_flagging: boolean, optional
            set this to use both pieter and christian flagging
        eclipse_flagging: boolean, optional
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
        ec_orbits = numpy.zeros(n_eff_orbits, dtype=int)
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

        time_method_pieter    = numpy.zeros(1024, dtype=numpy.byte)
        time_method_christian = numpy.zeros(1024, dtype=numpy.byte)

        #
        # compute "RTS ranking" of darklimb histogram, based on empirical criteria
        #
        #print dl_rtslevel.shape
        dl_rtslice = numpy.transpose(dl_rtslevel[1,:,0]) # checks if 2nd mode is active (i.e. multi-mode, i.e. rts)
        eclipse_histo_method = numpy.zeros(1024, dtype=numpy.byte)
        sum_ec_peaks =numpy.sum(ec_npeaks,axis=1)
        eclipse_histo_method[795:1024] = sum_ec_peaks[0:229] > 3
        dl_rts = dl_rtslice > 0
        msk = -1
        darklimb_histo_method = numpy.array(((dl_rts*msk) & numpy.transpose((1 + (jumpratios >= 1) * (peakratios < 2)))), dtype=numpy.byte)

        combined_arr = numpy.zeros(1024)

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
            #stds  = numpy.reshape(dl_stds[0,pix,:],?)
            #means = numpy.reshape(dl_means[0,pix,:],?)

            # remove trailing zeroes
            idx = numpy.where(means == 0)
            if idx[0].size > 1:
                stds  = stds[0:idx[0]]
                means = means[0:idx[0]]
            else:
                if idx[0] >= 0:
                    continue

            dl_med_std = numpy.median(stds)
            dl_diffs   = (means-numpy.roll(means,1))[1:]
            idx_step   = numpy.where(numpy.abs(dl_diffs) > 2*dl_med_std)
            idx_noise  = numpy.where(stds > 2*dl_med_std)
            time_method_pieter[pix] = (idx_step[0].size > 0) or (idx_noise[0].size > 0)

            #
            # find possible rts in time domain (christian's method): 
            # - extreme range 
            # (lower and upper bounds check not possible because not dark corrected!) 
            #

            minmean   = numpy.min(means)
            maxmean   = numpy.max(means)
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
    print("Pixelmask unit test:")
    p = PixelQuality()
    orbit = 10000
    print("initialised.")
    p.calculate(orbit)
    print("orbit %d computed" % orbit)
    p.write()
    print("Pixel mask data for orbit %d written to db." % orbit)

