# -*- coding: iso-8859-1 -*-
# module containing pixelmask class.
import ConfigParser
import numpy
import matplotlib.pyplot as plt
import h5py
import logging
import datetime

from read_statedark_module import sdmf_read_statedark
from darklimb_io_module import sdmf_read_rts_darklimb_record

# class that handles pixelmask computation. should be a 1:1 copy of the SDMF 3.1
# IDL code.
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
        
        cfg_file = open('default3.1.cfg')
        try:
            self.cfg = self.get_config(cfg_file)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+'default3.1.cfg!')
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

            self.errorlc_thres[channelindex]    = (self.cfg['dc_err_thres'])[i]
            self.residual_thres[channelindex]   = (self.cfg['res_thres'])[i]
            self.exposuretime_max[channelindex] = (self.cfg['dc_sat_time'])[i]
            self.sig_max[channelindex]          = (self.cfg['max_signal'])[i]

    # load configuration from file
    def get_config(self, config_file):
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
        return dict

    # perform the actual pixel mask computation
    def calculate(self, orbit, test=None):
        
        self.orbit = orbit
        darkids = self.cfg['deadflaggingstates']

        #
        # load dark fit
        #

        fdark = h5py.File(self.cfg['db_dir']+self.cfg['dark_fname'], 'r')
        gdark = fdark['DarkFit']
        darkcurrent_dset  = gdark["darkCurrent"]
        darkcurrenterror_dset  = gdark["darkCurrentError"]
        analogoffset_dset = gdark["analogOffset"]
        darkcurrent      = darkcurrent_dset[orbit-1,:]
        darkcurrenterror = darkcurrenterror_dset[orbit-1,:]
        analogoffset     = analogoffset_dset[orbit-1,:]

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
                logging.warning("oh noes, no dark data for orbit %d!" % orbit)
                return
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
        # correct limits in case of OCR43 darks, to keep pixelmasks similar to 
        # pre-OCR43 days.
        # different for every Epitaxx channel: ch6+ seems to deviate from ch8. 
        # ch7, which suffers from light leakage and for which the pixelmask has 
        # never been properly defined, is wildly different since OCR43. since 
        # nobody uses ch7, i'll try to adapt bad pixel count to pre-OCR43 
        # levels, but not meticulously. 
        #

        if orbit >= 43362:
            residual_thres[6*1024:7*1024-1]  *= .75
            exposuretime_max[6*1024:7*1024-1] *= 1.33
            errorlc_thres[6*1024:7*1024-1] *= .4
            # ch6+
            residual_thres[5*1024+795:6*1024-1] *= .2

        #
        # check invalid data points in the dark current calculation
        # if the values are 0, the pixel can not be corrected and is 
        # therefore unusable
        #

        invalid = (analogoffset == 0) | (darkcurrent == 0)

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
        #       maxsig: max signal per channel excluding AO and LC. over the 
        #               sahara, BU/s
        # TODO: this is wrong. you use the exposure time from the non-equator states, with signal values from the equator states.. 
        # you need to use signal values from the non-equator states. 
        # in practice this would mean 58000 BU/s limit instead of 55000 BU/s (i guess!)

        maxadc     = 65535
        lowlc      = - analogoffset / self.exposuretime_max
        highlc     = (maxadc-analogoffset)/self.exposuretime_max - self.sig_max
        darkcursat = (darkcurrent < lowlc) | (darkcurrent > highlc)

        #
        # Check for violation of the derived dark current error
        # The error in these parameters is not allowed to exceed a certain 
        # fraction of the parameter itself
        #

        #idx = numpy.where(darkcurrent > 0) # probably div-by-0 in here
        #print idx
        darkerrorfrac = abs(darkcurrenterror / darkcurrent) # TODO: store this
        dc_err_mask = darkerrorfrac > self.errorlc_thres

        #
        # check whether the darkcurrent residuals exceed their error
        #

        numdark = len(darkids)
        tmp     = numpy.zeros(self.pixelcount, dtype=bool)
        for darkid in darkids:
            state_string = "State_"+format(darkid, "02d")

            #
            # correct dark measurements with dark fit (sounds way too funny)
            #
            
            pets = readoutpet[state_string]
            corrmean  = readoutmean[state_string] 
            corrmean -= darkcurrent * pets + analogoffset
            corrnoise = readoutnoise[state_string]

            i_row = 0
            for meanrow in corrmean:
                noiserow = corrnoise[i_row]
                #phase = correcteddata[index[j]].phase
                #if phase > 0.0 and phase < 0.3:
                goodnoise=numpy.where(noiserow > 0)
                # if there are invalid noise figures, generate an error!
                if goodnoise[0].size == 0:
                    print 'invalid noise(=0) in residual criterion!'
                    return
                tmp = tmp | (abs(meanrow / noiserow) > self.residual_thres)
                i_row+=1
                
        residual = tmp

        #
        # compute RTS mask
        #
        
        rts = numpy.zeros(self.pixelcount)
        rts6 = numpy.zeros(self.numpixels)
        rts8 = numpy.zeros(self.numpixels)
        rts6 = self.calculate_rts_rank(orbit, channel=6, pieter_flagging=True, 
            christian_flagging=False, eclipse_flagging=True)
        rts8 = self.calculate_rts_rank(orbit, channel=8, pieter_flagging=True, 
            christian_flagging=False, eclipse_flagging=True)
        rts[5*1024:6*1024] = rts6
        rts[7*1024:8*1024] = rts8
        
        #
        # combine the masks using logical OR
        #

        combined = (rts > 0) | darkcursat | invalid | dc_err_mask | residual
        
        self.rts_mask = rts
        self.combined_mask = numpy.array(combined, dtype='=u1')
        self.darkcursat_mask = numpy.array(darkcursat, dtype='=u1')
        self.dc_err_mask = numpy.array(dc_err_mask, dtype='=u1')
        self.residual_mask = numpy.array(residual, dtype='=u1')
        self.invalid_mask = numpy.array(invalid, dtype='=u1')

        
    # creates a new mask array in the database
    def create_dbmaskarray(self, f, name, orbit):
        print "creating "+name
        dims = (self.maxorbits,self.pixelcount)
        if (orbit < 1) or (orbit > self.maxorbits):
            raise ValueError('orbit %d out of range [1,%d]' %
                             (orbit,self.maxorbits))
        dtype = numpy.dtype('=u1')
        dat = numpy.zeros(dims, dtype=dtype)
        f.create_dataset(name, dims, dtype=dtype, chunks=(16,self.pixelcount), 
                         compression='gzip', compression_opts=3, 
                         maxshape=[None,self.pixelcount], data=dat)

    # write data for this orbit to database
    def write(self):
        
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
            print "creating db "+db_fname+"..."
            f = h5py.File(db_fname,'w')
            print 'creating metaTable'
            metadata = numpy.empty((self.maxorbits), dtype=meta_dtype)
            f.create_dataset("metaTable", (self.maxorbits,), dtype=meta_dtype, 
                             chunks=(1024,), compression='gzip', 
                             compression_opts=3, maxshape=None, data=metadata)
            self.create_dbmaskarray(f, "RTS", orbit)
            self.create_dbmaskarray(f, "combined", orbit)
            self.create_dbmaskarray(f, "darkCurrentError", orbit)
            self.create_dbmaskarray(f, "darkCurrentSat", orbit)
            self.create_dbmaskarray(f, "residual", orbit)
            self.create_dbmaskarray(f, "invalid", orbit)
            f.close()
            print 'created db'

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
            logging.error("failed to open database %s" % db_fname)

    #
    # PURPOSE:
    #       generate pixelmasks for specified orbit list and channel
    #
    # INPUT:
    #       orbit_list: 1D array of orbits
    #
    # RETURNS:
    #       combined RTS mask
    #
    # KEYWORDS:
    #       channel: input, scalar, integer : 6 or 8, 6 implies channel 6+
    #        default = 6
    #       test: set this to plot bar graph instead of writing to db.
    #       pieter_flagging: set this to use time domain flagging, "pieter" style 
    #       christian_flagging: set this to use time domain flagging, "christian" 
    #        style 
    #       timedomain_flagging: set this to use both pieter and christian flagging
    #       eclipse_flagging: set this to use histogram eclipse data  
    #
    # NOTES: 
    #       - time domain flagging is very strict
    #       - eclipse flagging involves 3 orbits (!) because this is required to get
    #       decent statistics. this is probably a bit too strict. 
    #
    def calculate_rts_rank(self, orbit, channel=6, test=False, 
        pieter_flagging=True, christian_flagging=False, 
        timedomain_flagging=False, eclipse_flagging=False):

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
            print 'error calling sdmf_read_statedark_!'
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
            print 'got eclipse rts data'

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
            print 'error calling sdmf_read_rts_darklimb_record ['+str(status)+']!'
            return

        dl_rtslevel = data['rtslevel']
        dl_stds = data['execstddev']
        dl_means = data['execmean']
        peakratios = data['peakratios']
        jumpratios = data['jumpratios']

        if test:
            print 'got darklimb rts data'

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
    print "Pixelmask unit test:"
    p = Pixelmask()
    orbit = 14000
    print "initialised."
    p.calculate(orbit)
    print "orbit %d computed" % orbit
    print "good pixel fractions per criterion for channel 8, orbit %d:" % orbit
    print "combined", numpy.mean(p.combined_mask[7*1024:8*1024-1])
    print "darkcurrent saturation", \
          numpy.mean(p.darkcursat_mask[7*1024:8*1024-1])
    print "darkcurrent error", numpy.mean(p.dc_err_mask[7*1024:8*1024-1])
    print "residual", numpy.mean(p.residual_mask[7*1024:8*1024-1])
    print "invalid", numpy.mean(p.invalid_mask[7*1024:8*1024-1])
    print "RTS", numpy.mean(p.rts_mask[7*1024:8*1024-1] > 0)
    p.write()
    print "Pixel mask data for orbit %d written to db." % orbit
