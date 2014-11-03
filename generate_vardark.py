#!/usr/bin/env python
#
from __future__ import print_function, division

import numpy as np
import h5py

#-- functions ------------------------------------------------------------------

class VarDarkdb:
    '''
    '''
    def create(self, h5_name, sz_phase=50, sz_channel=1024):
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
                                        dtype="f", chunks=chunks )
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

        return

    def open(self, h5_name, sz_phase=50, sz_channel=1024):
        self.h5_name = np.string_(h5_name)

        # TODO: check if this matches the existing db's dimensions
        self.channel_sz = sz_channel
        self.phase_sz = sz_phase

        self.fid = h5py.File(h5_name, "r+")
        self.ds_orbit = self.fid["dim_orbit"]
        self.ds_vdark = self.fid["varDark"]
        self.ds_err_lcs = self.fid["errorLCs"]
        self.ds_err_trends = self.fid["errorTrends"]
        self.n_write = self.ds_orbit.size

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

    def replace(self, orbit, vdark, err_trends=None, err_lcs=None):
        '''
        Replace variable dark values in database
        '''
        indx = np.argwhere(self.ds_orbit[:] == orbit)
        print(self.ds_orbit[:], orbit, indx)
        if indx[0].size > 0:
            self.ds_vdark[indx,:,:] = vdark
            if err_trends is not None:
                self.ds_err_trends[indx,:] = err_trends
            if err_lcs is not None:
                self.ds_err_lcs[indx,:] = err_lcs

    def append(self, orbit, vdark, err_trends=None, err_lcs=None):
        '''
        Append new variable dark values in database
        '''
        self.ds_orbit.resize( (self.n_write+1,) )
        self.ds_orbit[self.n_write] = orbit
        self.ds_vdark.resize( self.n_write+1, axis=0 )
        self.ds_vdark[self.n_write,:,:] = vdark
        if err_trends is not None:
            self.ds_err_trends.resize( self.n_write+1, axis=0 )
            print(err_trends.shape)
            self.ds_err_trends[self.n_write,:] = err_trends
        if err_lcs is not None:
            self.ds_err_lcs.resize( self.n_write+1, axis=0 )
            self.ds_err_lcs[self.n_write,:] = err_lcs
        self.n_write += 1

    def close( self ):
        self.fid.close()

def parseOrbitList(str):
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

def generate_vardark(vddb, ad, input_dbname, first_orbit, last_orbit, pixnr=None, useTrendFit=True):
    '''
    generate vardark product for specified orbits and from specified input database.

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
        flag to set more precise fit of daily variation trending
    '''

    #
    # open interpolated monthlies
    #

    fin = h5py.File(input_dbname, "r")
    orbit_range = [first_orbit, last_orbit]
    in_orblist = fin["orbits"]
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
    # get the raw darks
    #

    print("get darks..")
    if first_orbit < 43362 and last_orbit >= 43362:
        print("lump upto 43361")
        ad.lump([first_orbit, 43361])
        print("lump from 43362")
        ad.lump([43362, last_orbit])
    else:
        print("lump")
        ad.lump(orbit_range)
    print("finalize")
    ad.finalize()
    print("get data")
    n_darks, dummy, pets, coadds, readouts, noise, ephases = ad.get_range(orbit_range, autoLump=False)
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

    print("extending borders..")
    ephasesi = ephases.astype(np.int32)
    print(ephasesi, xmin, xmax)
    idxl = ephasesi == int(xmin)
    n_l = np.sum(idxl)
    idxr = ephasesi == int(xmax)
    n_r = np.sum(idxr)
    #print(n_l, n_r)
    ephases = np.concatenate((ephases[idxl]-1., ephases, ephases[idxr]+1.))
    pets = np.concatenate((pets[idxl], pets, pets[idxr]))
    coadds = np.concatenate((coadds[idxl], coadds, coadds[idxr]))
    readouts = np.concatenate((readouts[idxl,:], readouts, readouts[idxr,:]))
    noise = np.concatenate((noise[idxl,:], noise, noise[idxr,:]))
    print("done.")

    #
    # subtract interpolated analog offset to get thermal signal, and normalize by time
    #

    print("get thermal background signal..")
    thermal_background = readouts
    i_orbit = 0
    m = 0
    eorbits = ephases.astype(np.int32)
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

    print("compute trending points..")
    i_trend = 0
    i_orbit = 0
    avg_phi = 0.
    xt = np.array([trending_phase]) # trending_phase
    lcs = np.ones(n_pix) * 5000 # initial guess for thermal background signal. 5000 BU/s is a good average
    orbrange = range(int(np.min(ephases)), int(np.max(ephases)))
    n_tpts = len(orbrange)
    trending_phis = np.empty(n_tpts)
    trending_ys = np.empty([n_tpts, n_pix])
    err_lcs = np.empty([n_tpts, n_pix])
    err_trends = np.empty([n_tpts, n_pix])

    for orbit in orbrange:
        #print(orbit)
        if useTrendFit:
            # use least-means fit to find lc and trend. slower, but more accurate. 
            # this works because all parameters are now fixed except lc and trend. 
            aos = inter_aos[i_orbit,:]
            amps = inter_amps[i_orbit,:]
            channel_phase1 = inter_phases[i_orbit,0]
            channel_phase2 = inter_phases[i_orbit,1]
            channel_amp2 = inter_amp2[i_orbit]

            #x__, lcs_, res_trends, dum1, dum2 = fit_eclipse_orbit(ad, orbit, aos, lcs, amps, channel_amp2, channel_phase1, channel_phase2)
            list_ = fit_eclipse_orbit(ad, orbit, aos, lcs, amps, channel_amp2, channel_phase1, channel_phase2, give_errors=True)
            x__, lcs_, res_trends, err_lcs_, err_trends_, dum1, dum2 = list_
            err_lcs[i_trend, :] = err_lcs_
            err_trends[i_trend, :] = err_trends_

            for i_pix in range(n_pix):
                p = aos[i_pix], lcs_[i_pix], amps[i_pix], res_trends[i_pix], channel_phase1, channel_amp2, channel_phase2
                trending_ys[i_trend, i_pix] = scia_dark_fun2m(p, xt)
            avg_phi += trending_phase
            trending_phis[i_trend] = trending_phase+orbit
            i_trend += 1
        else:
            # faster, but less accurate
            idx = (ephases >= orbit) & (ephases < (orbit+1))
            local_phi = ephases[idx]
            local_y = thermal_background[idx,:]
            abs_dist1 = np.abs(local_phi - (trending_phase+orbit))
            idx1 = (np.argsort(abs_dist1))[0:6]
            phi1 = np.mean(local_phi[idx1])
            avg1 = np.mean(local_y[idx1,:], axis=0)
            trending_phis[i_trend] = phi1
            trending_ys[i_trend,:] = avg1
            if not np.isnan(phi1):
                avg_phi += (phi1%1)
                i_trend += 1
        i_orbit += 1

    avg_phi /= i_trend
    print(avg_phi, i_trend)
    trending_phis = trending_phis[0:i_trend]
    trending_ys = trending_ys[0:i_trend]
    print("done.")

    print("remove invalid entries..")
    idx_goodphi = np.isfinite(trending_phis)
    trending_phis = trending_phis[idx_goodphi]
    trending_ys = trending_ys[idx_goodphi,:]
    # filter a single pixel. probably a bad idea
    #idx_goody = np.isfinite(trending_ys[:,pixnr])
    #trending_phis = trending_phis[idx_goody]
    #trending_ys = trending_ys[idx_goody,:]
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
    if pixnr is not None:
        f = interp1d(trending_phis, trending_ys[:,pixnr])
        f2 = interp1d(trending_phis, trending_ys[:,pixnr], kind='cubic')
    else:
        print(trending_phis.shape, trending_ys.shape)
        #print(trending_phis)
        f = interp1d(trending_phis, trending_ys, axis=0)
    print("done.")

    print("interpolate trend..")
    xt = np.array([phi_t]) # trending_phase
    wave = np.empty([pts_per_orbit, n_pix])
    full_wave = np.empty(xnew.size)
    out_orblist = np.array([], dtype=np.int32)

    i = 0
    for orbit in range(int(ximin), int(ximax)):
        out_orblist = np.append(out_orblist, orbit)
        print(orbit)
        i_orbit = (np.where(in_orblist[:] == orbit))[0][0]
        aos = inter_aos[i_orbit,:]
        amps = inter_amps[i_orbit,:]
        channel_phase1 = inter_phases[i_orbit,0]
        channel_phase2 = inter_phases[i_orbit,1]
        channel_amp2 = inter_amp2[i_orbit]
        idx = xnewi == orbit
        xnew_ = xnew[idx]
        if pixnr is not None:
            p = aos[pixnr], lcs[pixnr], amps[pixnr], 0, channel_phase1, channel_amp2, channel_phase2
            wave_ = scia_dark_fun2n(p, xnew_) - scia_dark_fun2n(p, xt)
            full_wave[idx] = wave_ + f(xnew_)
        else:
            for i_pix in range(n_pix):
                p = aos[i_pix], lcs[i_pix], amps[i_pix], 0, channel_phase1, channel_amp2, channel_phase2
                wave_ = scia_dark_fun2n(p, xnew_) - scia_dark_fun2n(p, xt)
                wave[pts_per_orbit-xnew_.size:pts_per_orbit, i_pix] = wave_ 
            wave[pts_per_orbit-xnew_.size:pts_per_orbit, :] += f(xnew_)
            if vddb.new_entry(orbit):
                if useTrendFit:
                    vddb.append(orbit, wave, err_trends=err_trends[i,:], err_lcs=err_lcs[i,:])
                else:
                    vddb.append(orbit, wave)
            else:
                if useTrendFit:
                    vddb.replace(orbit, wave, err_trends=err_trends[i,:], err_lcs=err_lcs[i,:])
                else:
                    vddb.append(orbit, wave)
        i += 1

    fin.close() # close interpolated monthlies databse, won't be needing it anymore
    print("done.")

    if pixnr is not None:
        # single pixel
        import matplotlib.pyplot as plt
        plt.ticklabel_format(useOffset=False)
        print(plot_x.shape, plot_y.shape)
        print(trending_phis.shape, trending_ys.shape)
        plt.plot(plot_x,plot_y,'v', trending_phis, trending_ys[:,pixnr], 'o',xnew,f(xnew),'-', xnew, f2(xnew),'--', xnew, full_wave,'-')
        plt.legend(['orig data', 'avg data', 'linear', 'cubic', 'reconstruct'], loc='best')
        plt.show()

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import re
    import argparse
    from argparse import ArgumentParser, ArgumentTypeError
    import subprocess
    import matplotlib.pyplot as plt
    import warnings
    from scipy.interpolate import interp1d
    from scia_dark_functions import scia_dark_fun2n, scia_dark_fun2m
    #warnings.simplefilter("error") # warnings to errors
    #np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    from sciamachy_module import get_darkstateid, petcorr, n_chanpix
    from vardark_module import AllDarks, trending_phase, fit_eclipse_orbit
    import envisat

    n_pix = n_chanpix
    pts_per_orbit = 50
    #path = "/array/slot0B/SDMF/3.1/pieter/"
    path = ""
    dbname_long = path+"vardark_long.h5"
    dbname_short = path+"vardark_short.h5"

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

    if args.shortMode:
        print("PETs [0.0625, 0.125]")
        ad = AllDarks([0.0625, 0.125])
        dbname = dbname_short
        input_dbname = path+"interpolated_monthlies_short.h5"
    else:
        print("PETs [1.0, 0.5]")
        ad = AllDarks([1.0, 0.5])
        dbname = dbname_long
        input_dbname = path+"interpolated_monthlies_long.h5"

    if args.output_fname:
        dbname = args.output_fname
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
        first_orbit = 5000
        last_orbit = 55000
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
            interpolate_monthlies(input_dbname, path+"monthly_fits_short.h5")
        else:
            interpolate_monthlies(input_dbname, path+"monthly_fits_long.h5")

    print(input_dbname, "-> (generate_vardark) -> ", dbname)

    #
    # open or create output database 
    #

    vddb = VarDarkdb(dbname, sz_phase=pts_per_orbit, sz_channel=n_pix)

    #
    # do the work
    #

    generate_vardark(vddb, ad, input_dbname, first_orbit, last_orbit, pixnr=pixnr, useTrendFit=True)

    #
    # close down the output database
    #

    print("closing..")
    vddb.close()    
    print("done.")
