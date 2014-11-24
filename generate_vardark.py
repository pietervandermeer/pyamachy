#!/usr/bin/env python
#
from __future__ import print_function, division

import numpy as np
import h5py

#-- functions ------------------------------------------------------------------

class VarDarkdb:
    '''
    '''
    def create(self, h5_name, sz_phase=50, sz_channel=1024, use_trendfit=False, short_pet=False):
        '''
        initialize HDF5 file
        '''
        self.h5_name = np.string_(h5_name)
        self.n_write = 0
        self.channel_sz = sz_channel
        self.phase_sz = sz_phase

        # set cache-size and version bounding
        propfaid = h5py.h5p.create( h5py.h5p.FILE_ACCESS )
        settings = list( propfaid.get_cache() )
        settings[2] = 4 * 1024 * 1024          ## = 4 MB
        propfaid.set_cache( *settings )

        propfaid.set_libver_bounds( h5py.h5f.LIBVER_LATEST,
                                    h5py.h5f.LIBVER_LATEST )
        
        # create HDF5 file
        self.fp = h5py.h5f.create( self.h5_name, fapl=propfaid )
        self.fid = h5py.File( self.fp )

        self.fid.attrs['trendfit'] = use_trendfit
        self.fid.attrs['short_pet'] = short_pet

        #
        # define dimensions and datasets in HDF5 file
        #
        pixels = np.arange( self.channel_sz, dtype='u2' )
        dset = self.fid.create_dataset( 'dim_pixel', data=pixels )
        dset.dims.create_scale( self.fid['dim_pixel'], 
                      'This is a netCDF dimension but not a netCDF variable.' )
        self.ds_pixel = dset

        phases = np.arange( self.phase_sz, dtype='f' ) / float(self.phase_sz)
        dset = self.fid.create_dataset( 'dim_phase', data=phases )
        dset.dims.create_scale( self.fid['dim_phase'], 
                      'This is a netCDF dimension but not a netCDF variable.' )
        self.ds_phase = dset

        dset = self.fid.create_dataset( "dim_orbit", (0,), maxshape=(None,),
                                        dtype="u2" )
        dset.dims.create_scale( self.fid['dim_orbit'], 
                      'This is a netCDF dimension but not a netCDF variable.' )
        self.ds_orbit = dset

        chunks = (1, self.ds_phase.size, self.ds_pixel.size,)
        creshape = (0, self.ds_phase.size, self.ds_pixel.size)
        maxshape = (None, self.ds_phase.size, self.ds_pixel.size,)
        dset = self.fid.create_dataset( "varDark", creshape,
                                        maxshape=maxshape, 
                                        dtype="f") #, chunks=chunks ) # chunking disabled as current libraries give deprecation error
        dset.dims[0].attach_scale( self.fid['dim_orbit'] )
        dset.dims[0].label = 'absolute orbit number'
        dset.dims[1].attach_scale( self.fid['dim_phase'] )
        dset.dims[1].label = 'orbit phase (eclipse)'
        dset.dims[2].attach_scale( self.fid['dim_pixel'] )
        dset.dims[2].label = 'pixel ID'
        dset.attrs['long_name'] = np.string_("variable dark (channel 8)")
        dset.attrs['units'] = np.string_("BU/s")
        dset.attrs['description'] = \
            np.string_("SCIA variable dark signal of channel 8")
        self.ds_vdark = dset

        self.ds_err_lcs = self.fid.create_dataset("errorLCs", (0,sz_channel), maxshape=(None,sz_channel), dtype=np.float64)
        self.ds_err_trends = self.fid.create_dataset("errorTrends", (0,sz_channel), maxshape=(None,sz_channel), dtype=np.float64)
        self.ds_datapoint_counts = self.fid.create_dataset("dataPointCounts", (0,), maxshape=(None,), dtype=np.int16)
        self.ds_uncertainties = self.fid.create_dataset("uncertainties", (0,sz_channel), maxshape=(None,sz_channel), dtype=np.float64)

        return

    def open(self, h5_name, sz_phase=50, sz_channel=1024, use_trendfit=False, short_pet=False):
        self.h5_name = np.string_(h5_name)

        # TODO: check if this matches the existing db's dimensions
        self.channel_sz = sz_channel
        self.phase_sz = sz_phase

        try:
            self.fid = h5py.File(h5_name, "r+")
        except IOError:
            raise Exception("unable to open hdf5 file "+h5_name)

        if self.fid.attrs['trendfit'] != use_trendfit:
            raise Exception("hdf5 file "+h5_name+" has wrong attributes.")
        if self.fid.attrs['short_pet'] != short_pet:
            raise Exception("hdf5 file "+h5_name+" has wrong attributes.")

        try:
            self.ds_orbit = self.fid["dim_orbit"]
            self.ds_pixel = self.fid["dim_pixel"]
            self.ds_phase = self.fid["dim_phase"]
            self.ds_vdark = self.fid["varDark"]
            self.ds_err_lcs = self.fid["errorLCs"]
            self.ds_err_trends = self.fid["errorTrends"]
            self.ds_datapoint_counts = self.fid["dataPointCounts"]
            self.ds_uncertainties = self.fid["uncertainties"]
        except KeyError:
            raise Exception("unable to open one or more datasets in "+h5_name)

        self.n_write = self.ds_orbit.size
        n_pixel = self.ds_pixel.size
        n_phase = self.ds_phase.size
        if n_pixel != sz_channel:
            raise Exception("nr of pixels does not match between user-specified and existing database!")
        if n_phase != sz_phase:
            raise Exception("nr of phases does not match between user-specified and existing database!")

        return

    def __init__( self, h5_name, **kwargs):

        #
        # create if database doesn't exist 
        #

        if not h5py.is_hdf5(h5_name):
            self.create(h5_name, **kwargs)
            return

        #
        # otherwise just open existing
        #

        self.open(h5_name, **kwargs)

        return

    def new_entry( self, orbit ):
        '''
        Check if data is available for given orbit
        '''
        indx = np.argwhere( self.ds_orbit[:] == orbit )
        if indx.size > 0:
            return False
        else:
            return True

    def replace(self, orbit, vdark, err_trends=None, err_lcs=None, datapoint_count=None, uncertainties=None):
        '''
        Replace variable dark values in database
        '''
        boolindex = self.ds_orbit[:] == orbit
        indx = np.argwhere(boolindex)
        if indx[0].size == 0:
            return

        self.ds_vdark[indx,:,:] = vdark
        if err_trends is not None:
            self.ds_err_trends[indx,:] = err_trends
        if err_lcs is not None:
            self.ds_err_lcs[indx,:] = err_lcs
        if datapoint_count is not None:
            self.ds_datapoint_counts[indx[0][0]] = datapoint_count
        if uncertainties is not None:
            self.ds_uncertainties[indx,:] = uncertainties

    def append(self, orbit, vdark, err_trends=None, err_lcs=None, datapoint_count=None, uncertainties=None):
        '''
        Append new variable dark values in database
        '''
        self.ds_orbit.resize( (self.n_write+1,) )
        self.ds_orbit[self.n_write] = orbit
        self.ds_vdark.resize( self.n_write+1, axis=0 )
        self.ds_vdark[self.n_write,:,:] = vdark
        if err_trends is not None:
            self.ds_err_trends.resize( self.n_write+1, axis=0 )
            self.ds_err_trends[self.n_write,:] = err_trends
        if err_lcs is not None:
            self.ds_err_lcs.resize( self.n_write+1, axis=0 )
            self.ds_err_lcs[self.n_write,:] = err_lcs
        if datapoint_count is not None:
            self.ds_datapoint_counts.resize((self.n_write+1,))
            self.ds_datapoint_counts[self.n_write] = datapoint_count
        if uncertainties is not None:
            self.ds_uncertainties.resize( self.n_write+1, axis=0 )
            self.ds_uncertainties[self.n_write,:] = uncertainties
        self.n_write += 1

    def store(self, orbit, vdark, **kwargs):
        '''
        Store new variable dark values in database. Automatically checks whether to replace or append.
        '''
        if self.new_entry(orbit):
            self.append(orbit, vdark, **kwargs)
        else:
            self.replace(orbit, vdark, **kwargs)
        return

    def close( self ):
        self.fid.close()

def generate_vardark(vddb, ad, input_dbname, first_orbit, last_orbit, pixnr=None, useTrendFit=True):
    """
    Generate vardark product for specified orbits and from specified input database.

    Parameters
    ----------

    vddb : VarDarkdb object
        used to store the data as hdf5
    ad : AllDarks object 
        contains dark data for specified orbit range
    input_dbname : string
        name of interpolated monthlies database
    first_orbit : int 
        first orbit in orbit range
    last_orbit] : int 
        last orbit in orbit range
    pixnr : int, optional 
        [0..1023]. if not None will plot a pixel of choice instead of output data to database
    useTrendFit : boolean, optional
        flag that replaces simple trending procedure by fit (more accurate, but also more costly)
    """

    import logging

    #
    # open interpolated monthlies
    #

    fin = h5py.File(input_dbname, "r")
    in_orblist = fin["orbits"]
    n_orbits_in = in_orblist.size
    inter_aos = fin["aos"]
    inter_amps = fin["amps"]
    inter_phases = fin["phases"]
    inter_amp2 = fin["amp2"]

    # plot phaseshift vs orbit over the entire mission
    # import matplotlib.pyplot as plt
    # plt.ticklabel_format(useOffset=False)
    # plt.plot(in_orblist[:],inter_phases[:,0])
    # plt.show()

    #
    # handle orbit range
    #

    first_calib_orbit = min(in_orblist[:])
    if first_orbit is None:
        first_orbit = first_calib_orbit

    last_calib_orbit = max(in_orblist[:])
    if last_orbit is None:
        last_orbit = last_calib_orbit

    orbit_range = [first_orbit, last_orbit]

    if first_calib_orbit > first_orbit:
        raise Exception("first orbit before first montly calibration orbit "+str(first_calib_orbit))
    if last_calib_orbit < last_orbit:
        raise Exception("last orbit after last montly calibration orbit "+str(last_calib_orbit))

    print("orbit range =", orbit_range)

    #
    # get the raw darks
    #

    print("get darks..")
    n_darks, dummy, pets, coadds, readouts, noise, ephases = ad.get_range(orbit_range)
    eorbits = ephases.astype(np.int32)
    print("done.")

    #
    # sort darks
    #

    print("sort darks..")
    idx = np.argsort(ephases)
    ephases = ephases[idx]
    pets = pets[idx]
    coadds = coadds[idx]
    readouts = readouts[idx,:]
    noise = noise[idx,:]
    xmin = np.min(ephases)
    xmax = np.max(ephases)
    print("done.")

# may be useless.. as trending limits itself to single points per orbit
#    print("extending borders..")
#    ephasesi = ephases.astype(np.int32)
#    print(ephasesi, xmin, xmax)
#    idxl = ephasesi == int(xmin)
#    n_l = np.sum(idxl)
#    idxr = ephasesi == int(xmax)
#    n_r = np.sum(idxr)
#    ephases = np.concatenate((ephases[idxl]-1., ephases, ephases[idxr]+1.))
#    pets = np.concatenate((pets[idxl], pets, pets[idxr]))
#    coadds = np.concatenate((coadds[idxl], coadds, coadds[idxr]))
#    readouts = np.concatenate((readouts[idxl,:], readouts, readouts[idxr,:]))
#    noise = np.concatenate((noise[idxl,:], noise, noise[idxr,:]))
#    print("done.")

    #
    # subtract interpolated analog offset to get thermal signal, and normalize by time
    #

    print("get thermal background signal..")
    thermal_background = readouts
    i_orbit = 0
    m = 0
    for orbit in in_orblist[:]:
        aos = inter_aos[i_orbit,:] 
        n = np.sum(eorbits == orbit)
        if n > 0:
            thermal_background[m:m+n,:] -= (np.matrix(aos).T * np.ones(n)).T
        i_orbit += 1
        m += n
    thermal_background /= np.matrix(pets).T * np.ones(n_pix)
    print("done.")

    if pixnr is not None:
        plot_x = ephases
        plot_y = thermal_background[:,pixnr]

    #
    # determine trending point for each orbit (including the extended borders)
    #

    print("compute trends and lc offsets..")
    i_trend = 0
    i_orbit = 0
    avg_phi = 0.
    xt = np.array([trending_phase]) # trending_phase
    lcs = np.ones(n_pix) * 5000 # initial guess for thermal background signal. 5000 BU/s is a good average
    orbrange = range(int(np.min(ephases)), int(np.max(ephases)))
    n_tpts = len(orbrange)
    # TODO: put phi, orbit, y, err_lc, err_trend, uncertainties into dict.. otherwise stoo much of a hassle when sorting or filteringm
    trending_phis = np.empty(n_tpts)
    trending_orbits = np.empty(n_tpts)
    trending_ys = np.empty([n_tpts, n_pix])
    err_lcs = np.empty([n_tpts, n_pix])
    err_trends = np.empty([n_tpts, n_pix])
    datapoint_count = np.empty([n_orbits_in]) # may be used for quality estimation (2 or less datapoints will make a fit impossible, even)
    uncertainties = np.empty([n_tpts, n_pix])

    for orbit in orbrange:
        print(orbit)
        if useTrendFit:
            # use least-means fit to find lc and trend. slower, but more accurate. 
            # this works because all parameters are now fixed except lc and trend. 
            idx = in_orblist[:] == orbit
            if np.sum(idx) > 0:
                idx = np.where(idx)[0][0]
                aos = inter_aos[idx,:]
                amps = inter_amps[idx,:]
                channel_phase1 = inter_phases[idx,0]
                channel_phase2 = inter_phases[idx,1]
                channel_amp2 = inter_amp2[idx]

                try:
                    #x__, lcs_, res_trends, dum1, dum2 = fit_eclipse_orbit(ad, orbit, aos, lcs, amps, channel_amp2, channel_phase1, channel_phase2)
                    list_ = fit_eclipse_orbit(ad, orbit, aos, lcs, amps, channel_amp2, channel_phase1, channel_phase2, give_errors=True, verbose=False)
                    x__, lcs_, res_trends, err_lcs_, err_trends_, dum1, dum2, uncertainty = list_
                    err_lcs[i_trend, :] = err_lcs_
                    err_trends[i_trend, :] = err_trends_

                    datapoint_count[i_orbit] = x__[0].size
                    uncertainties[i_trend,:] = uncertainty
                    for i_pix in range(n_pix):
                        p = aos[i_pix], lcs_[i_pix], amps[i_pix], res_trends[i_pix], channel_phase1, channel_amp2, channel_phase2
                        trending_ys[i_trend, i_pix] = scia_dark_fun2m(p, xt)
                    avg_phi += trending_phase
                    trending_phis[i_trend] = trending_phase+orbit
                    trending_orbits[i_trend] = orbit
                    i_trend += 1
                except Exception as e:
                    # just skip to next orbit and don't store any data, but do log a warning!
                    logging.warning("failed to fit orbit.. "+str(e))
                    datapoint_count[i_orbit] = 0         
        else:
            # faster, but less accurate
            idx = (ephases >= orbit) & (ephases < (orbit+1))
            datapoint_count[i_orbit] = idx.size
            local_phi = ephases[idx]
            local_y = thermal_background[idx,:]
            print("local_y=", local_y[:,pixnr])
            abs_dist1 = np.abs(local_phi - (trending_phase+orbit))
            idx1 = (np.argsort(abs_dist1))[0:6]
            phi1 = np.mean(local_phi[idx1])
            avg1 = np.mean(local_y[idx1,:], axis=0)
            print("avg1=", avg1[pixnr])
            trending_phis[i_trend] = phi1
            trending_ys[i_trend,:] = avg1
            trending_orbits[i_trend] = orbit
            datapoint_count[i_trend] = idx1.size #TODO: test.. idx1 isn't a tuple of np arrays?s
            if not np.isnan(phi1):
                avg_phi += (phi1%1)
                i_trend += 1
        i_orbit += 1

    avg_phi /= i_trend
    print(avg_phi, i_trend)
    trending_phis = trending_phis[0:i_trend]
    trending_ys = trending_ys[0:i_trend]
    uncertainties = uncertainties[0:i_trend,:]
    trending_orbits = trending_orbits[0:i_trend]
    err_lcs = err_lcs[0:i_trend,:]
    err_trends = err_trends[0:i_trend,:]
    datapoint_count = datapoint_count[0:i_trend]
    print("done.")

    print("remove invalid entries..")
    idx_goodphi = np.isfinite(trending_phis)
    trending_phis = trending_phis[idx_goodphi]
    trending_ys = trending_ys[idx_goodphi,:]
    uncertainties = uncertainties[idx_goodphi,:]
    trending_orbits = trending_orbits[idx_goodphi]
    err_lcs = err_lcs[idx_goodphi,:]
    err_trends = err_trends[idx_goodphi,:]
    datapoint_count = datapoint_count[idx_goodphi]
    # filter a single pixel. probably a bad idea
    #idx_goody = np.isfinite(trending_ys[:,pixnr])
    #trending_phis = trending_phis[idx_goody]
    #trending_ys = trending_ys[idx_goody,:]
    #...
    print("done.")

    ximin = int(xmin)
    ximax = int(xmax)
    xnew = np.linspace(ximin, ximax, (ximax-ximin)*pts_per_orbit+1)
    phi_t = avg_phi # trending_phase
    print("phi_t=",phi_t)

    #
    # interpolate.. then visualize a pixel of choice or store data
    #

    xnewi = xnew.astype(np.int32)
    print("generating interpolators..")
    # extend range for interpolation
    trending_phis = np.concatenate((trending_phis[:1].astype('i'),trending_phis,trending_phis[-1:].astype('i')+1))
    if pixnr is not None:
        trending_ys = np.concatenate((trending_ys[:1],trending_ys,trending_ys[-1:]))
        f = interp1d(trending_phis, trending_ys[:,pixnr])
        f2 = interp1d(trending_phis, trending_ys[:,pixnr], kind='cubic')
    else:
        trending_ys = np.concatenate((trending_ys[:1,:],trending_ys,trending_ys[-1:,:]))
        #print(trending_phis.shape, trending_ys.shape)
        print(trending_phis)
        f = interp1d(trending_phis, trending_ys, axis=0)
    print("done.")

    print("interpolate trend..")
    xt = np.array([phi_t]) # trending_phase
    wave = np.empty([pts_per_orbit, n_pix])
    full_wave = np.empty(xnew.size)
    out_orblist = np.array([], dtype=np.int32)

    i = 0
    dummy_uncertainties = np.empty([n_pix])
    dummy_uncertainties[:] = np.nan
    for orbit in range(int(ximin), int(ximax)):
        out_orblist = np.append(out_orblist, orbit)
        print(orbit)
        i_orbit = (np.where(in_orblist[:] == orbit))[0][0]
        res = (np.where(trending_orbits == orbit))[0]
        if res.size > 0:
            i_trend = res[0]
        else:
            i_trend = -1
        aos = inter_aos[i_orbit,:]
        amps = inter_amps[i_orbit,:]
        channel_phase1 = inter_phases[i_orbit,0]
        channel_phase2 = inter_phases[i_orbit,1]
        channel_amp2 = inter_amp2[i_orbit]
        idx = xnewi == orbit
        xnew_ = xnew[idx]
        if pixnr is not None:
            p = aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phase1, channel_amp2, channel_phase2
            wave_ = scia_dark_fun2n(p, xnew_) - scia_dark_fun2n(p, xt) # orbital variation wave only, no lc offset
            full_wave[idx] = wave_ + f(xnew_) # add interpolated lc offset (daily+seasonal variation)
        else:
            for i_pix in range(n_pix):
                p = aos[i_pix], lcs[i_pix], amps[i_pix], 0, channel_phase1, channel_amp2, channel_phase2
                wave_ = scia_dark_fun2n(p, xnew_) - scia_dark_fun2n(p, xt) # orbital variation wave only, no lc offset
                wave[pts_per_orbit-xnew_.size:pts_per_orbit, i_pix] = wave_ # add interpolated lc offset (daily+seasonal variation)
            #print(xnew_)
            wave[pts_per_orbit-xnew_.size:pts_per_orbit, :] += f(xnew_)
            if useTrendFit:
                if i_trend >= 0:
                    vddb.store(orbit, wave, datapoint_count=datapoint_count[i_trend], uncertainties=uncertainties[i_trend,:])
                else:
                    vddb.store(orbit, wave, datapoint_count=0, uncertainties=dummy_uncertainties)
            else:
                vddb.store(orbit, wave)
        i += 1

    fin.close() # close interpolated monthlies database, won't be needing it anymore
    print("done.")

    if pixnr is not None:
        # single pixel
        import matplotlib.pyplot as plt
        plt.ticklabel_format(useOffset=False)
        print(plot_x.shape, plot_y.shape)
        print(trending_phis.shape, trending_ys.shape)
        plt.title("Variable dark current for pixel "+str(pixnr)+", ch8")
        plt.xlabel("Orbit number")
        plt.ylabel("Dark current (BU/s)")
        plt.plot(plot_x,plot_y,'v', trending_phis, trending_ys[:,pixnr], 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--', xnew, full_wave,'-')
        plt.legend(['orig data', 'avg data', 'linear', 'cubic', 'reconstruct'], loc='best')
        plt.show()

    return

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser, ArgumentTypeError
    import subprocess
    import matplotlib.pyplot as plt
    import warnings
    from scipy.interpolate import interp1d
    from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
    import os.path
    #warnings.simplefilter("error") # warnings to errors
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    from sciamachy_module import get_darkstateid, petcorr, n_chanpix
    from vardark_module import AllDarks, trending_phase, fit_eclipse_orbit
    import envisat
    from envisat import parseOrbitList

    n_pix = n_chanpix
    pts_per_orbit = 50
    #path = "/array/slot0B/SDMF/3.1/pieter" # risky.. might lead to overwriting production database
    #path = "/SCIA/SDMF31/pieter"
    path = "./" # default.. safe 

    #
    # parse command line arguments
    #

    parser = argparse.ArgumentParser(description=
             'Computes channel 8 thermal background signal in BU/s for specified orbit range. Uses interpolated_monthlies.h5.')
    parser.add_argument('-o', '--output', dest='output_fname', type=str)
    parser.add_argument('--config', dest='config_file', type=file, 
                        default='default3.1.cfg')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true')
    parser.add_argument('-V', '--version', action='version', 
                        version='%(prog)s 0.1')
    parser.add_argument('-P', '--path', dest='path')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--orbitrange', default='all', 
                        help='sets orbit range f.e. "43000-44000", "all"', type=parseOrbitList)
    # parameters used especially for the plot
    parser.add_argument('-s', '--short', action='store_true', 
                        dest='shortMode', help="compute product for short PET's instead of long ones.")
    parser.add_argument('-p', '--pixnr', action='store', type=int, default=None,
                        dest='pixnr', help="pixel number to be examined [0..1023]")
    parser.add_argument('-f', '--fittrend', action='store_true', default=False,
                        dest='useFitTrend', help="least-squares fit instead of just getting average trending point.")
    args = parser.parse_args()

    #
    # handle command line arguments
    #

    if args.path is not None:
        path = args.path

    # override default output db name and path when output db specified!
    if args.output_fname:
        dbname = args.output_fname
        path = os.path.dirname(dbname)
        if not path:
            path = "."
    print("path=", path)

    if args.shortMode:
        print("PETs [0.0625, 0.125]")
        ad = AllDarks([0.0625, 0.125])
        dbname = path+"/vardark_short.h5"
        input_dbname = path+"/interpolated_monthlies_short.h5"
    else:
        print("PETs [1.0, 0.5]")
        ad = AllDarks([1.0, 0.5])
        dbname = path+"/vardark_long.h5"
        input_dbname = path+"/interpolated_monthlies_long.h5"
    print("dbname=", dbname)

    if args.pixnr is not None:
        if args.pixnr >= 0 and args.pixnr < n_pix:
            pixnr = args.pixnr
            print("pixel nr set, enabling plot mode (not generating database!).")
        else:
            print(args.pixnr, "out of range [0.."+str(n_pix-1)+"]!")
            quit()
        print("pixnr=", pixnr)
    pixnr = args.pixnr

    if args.orbitrange is None:
        first_orbit = None
        last_orbit = None
    else:
        first_orbit = args.orbitrange[0]
        last_orbit = args.orbitrange[1]
    print("orbitrange=", first_orbit, "-", last_orbit)

    # use fitting procedure for every orbit y/n. if no, just find a good trending point by averaging local darks (much faster).
    useTrendFit = args.useFitTrend
    print("useTrendFit =", useTrendFit)

    if not h5py.is_hdf5(input_dbname):
        from interpolate_monthlies import interpolate_monthlies
        print(input_dbname, "not present. generating...")
        if args.shortMode:
            interpolate_monthlies(input_dbname, path+"/monthly_fits_short.h5")
        else:
            interpolate_monthlies(input_dbname, path+"/monthly_fits_long.h5")

    print(input_dbname, "-> (generate_vardark) -> ", dbname)

    #
    # open or create output database 
    #

    if pixnr is not None:
        vddb = None
    else:
        vddb = VarDarkdb(dbname, 
                         sz_phase=pts_per_orbit, sz_channel=n_pix, 
                         use_trendfit=useTrendFit, short_pet=args.shortMode)

    #
    # do the work
    #

    generate_vardark(vddb, ad, input_dbname, first_orbit, last_orbit, pixnr=pixnr, useTrendFit=useTrendFit)

    #
    # close down the output database
    #

    print("closing..")
    vddb.close()    
    print("done.")
