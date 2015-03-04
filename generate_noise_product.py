#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
generate a simple noise product from sdmf readout noise.
the use of this generator is the fact that level 0 data are not always sliced equally for each orbit. 
instead of containing 15 dark executions they may contain 10 and the next 20.
this leads to variation from orbit to orbit.
this generator overcomes this by taking all data belonging to the right orbit (using eclipse definition).
the product also takes into account the changed dark definition (OCR-43) @ orbit 43362:
it will always have the right noise figure for the Pixel Exposure Time (PET) of choice.
"""

from __future__ import print_function, division

#import cPickle as pickle # python 2 speed-up
import pickle
import numpy as np 
import h5py

from sciamachy import get_darkstateid, read_extracted_states_ch8
from envisat import PhaseConverter

#-- globals --------------------------------------------------------------------

n_pix = 1024
use_sdmf_30 = True # use SDMF3.0 data? use SDMF3.1 otherwise
# storage method. interesting.. pickle is incompatible with ctypes.. so never mind.. 
#storage_method = 's' # pickle string
#storage_method = 'p' # pickle binary file
storage_method = 'h' # hdf5

#-- functions ------------------------------------------------------------------

class Noise:
    def __init__(self, pet):
        self.pet = pet
        if use_sdmf_30:
            self.calib_db = '/SCIA/SDMF30/sdmf_extract.h5' # SDMF 3.0
        else:
            self.calib_db = '/SCIA/SDMF31/sdmf_extract_calib.h5' # SDMF 3.1
        self.pc = PhaseConverter()

    # extract noise data from sdmf_extract_calib database
    def extract(self, orbit_range):

        #
        # read data
        #

        if orbit_range[0] < 43362 and orbit_range[1] >= 43362:
            pre_ocr43_sid = get_darkstateid(self.pet, 43361)
            post_ocr43_sid = get_darkstateid(self.pet, 43362)
            pre_range = [orbit_range[0], 43361]
            post_range = [43362, orbit_range[1]]
            d1 = read_extracted_states_ch8(pre_range, pre_ocr43_sid, self.calib_db, readoutNoise=True, errorInTheMean=False, readoutCount=True, columnOrder30=True)
            d2 = read_extracted_states_ch8(post_range, post_ocr43_sid, self.calib_db, readoutNoise=True, errorInTheMean=False, readoutCount=True, columnOrder30=True)
            c1 = d1['readoutCount'][:]
            c2 = d2['readoutCount'][:]
            self.count = np.concatenate((c1,c2))
            j1 = d1['mtbl'][:]['julianDay']
            j2 = d2['mtbl'][:]['julianDay']
            self.mjds = np.concatenate((j1,j2))
            n1 = d1['readoutNoise'][:,:]
            n2 = d2['readoutNoise'][:,:]
            self.noise = np.concatenate((n1,n2))
        else:
            stateid = get_darkstateid(self.pet, orbit_range[0])
            d = read_extracted_states_ch8(orbit_range, stateid, self.calib_db, readoutNoise=True, errorInTheMean=False, readoutCount=True, columnOrder30=True)
            self.count = d['readoutCount'][:]
            self.mjds = d['mtbl'][:]['julianDay']
            self.noise = d['readoutNoise'][:,:]

        #
        # deduplicate if necessary (shouldn't be)
        #

        # self.mjds, idx = np.unique(self.mjds, return_index=True)
        # print(idx)
        # self.count = self.count[idx]
        # self.noise = self.noise[idx]

        #
        # convert mjd to orbit to get reliable orbit numbers based on eclipse phase
        #

        ephases, self.orbits = self.pc.get_phase(self.mjds, eclipseMode=True, getOrbits=True)
        self.orbits += ephases # phase can be outside of range [0..1], so add to orbit and truncate back to int to get the real orbit nrs!
        self.orbits.view('int32')

        return

    # averages noise figures so that they are per orbit
    def finalize(self):
        # sort.. to group state executions by orbit
        idx = np.argsort(self.mjds)
        self.mjds = self.mjds[idx]
        self.orbits = self.orbits[idx]
        print(self.noise.shape)
        self.noise = self.noise[idx,:]
        self.count = self.count[idx,:]

        #
        # average the noise figure for each orbit
        #

        # get position and group size
        firsts = np.nonzero(np.diff(self.orbits))[0] + 1
        firsts = np.insert(firsts, 0, 0)
        diffs = np.diff(firsts)
        exec_counts = np.insert(diffs, -1, self.orbits.size - firsts[-1])

        # create output array (where all orbits are collapsed into single columns)
        n_orbits = exec_counts.size
        noise_out = np.empty([n_orbits, n_pix])
        self.meas_count = np.empty([n_orbits, n_pix])

        # average each orbit
        for i in range(n_orbits):
            count = exec_counts[i]
            first = firsts[i]
            sdev = self.noise[first:first+count, :]  # stddev
            meas_count = np.sum(self.count[first:first+count], axis=0)  # total count of all measurements for this orbit
            self.meas_count[i,:] = meas_count # store total measurement count for this orbit
            #noise_out[i,:] = np.sqrt(np.mean(np.square(sdev), axis=0) / meas_count) # root(mean_square / count) -> error in the mean
            noise_out[i,:] = np.sqrt(np.mean(np.square(sdev), axis=0)) # average stddev (RMS)

        self.noise = noise_out

        # we want only one entry per orbit in the database..
        self.orbits = np.unique(self.orbits)

        return

    def get_groupname(self):
        return "pet"+str(self.pet)

    def dump_hdf5(self, fname):
        f = h5py.File(fname, "a")
        grpn = self.get_groupname()
        n_dset = f.create_dataset(grpn+"/noise", self.noise.shape, dtype='f')
        n_dset.attrs["long_name"] = np.string_("Noise figures per orbit per pixel, for one pixel exposure time.")
        n_dset.attrs["units"] = np.string_("BU")
        n_dset.attrs["description"] = np.string_("Noise figures")
        n_dset[:,:] = self.noise
        c_dset = f.create_dataset(grpn+"/count", self.meas_count.shape, dtype='i')
        c_dset.attrs["long_name"] = np.string_("Number of measurements used for noise figures (per orbit per pixel).")
        c_dset.attrs["units"] = np.string_("-")
        c_dset.attrs["description"] = np.string_("Number of measurements used.")
        c_dset[:,:] = self.meas_count
        o_dset = f.create_dataset(grpn+"/orbits", self.orbits.shape, dtype='i')
        o_dset[:] = self.orbits
        o_dset.attrs["long_name"] = np.string_("Absolute orbit numbers.")
        o_dset.attrs["units"] = np.string_("-")
        o_dset.attrs["description"] = np.string_("Absolute orbit numbers.")
        f.close()

    def load_hdf5(self, fname):
        f = h5py.File(fname, "r")
        grpn = self.get_groupname()
        self.noise = f[grpn+"/noise"][:,:]
        self.orbits = f[grpn+"/orbits"][:]
        f.close()

#-- main -----------------------------------------------------------------------

def test_load_store():
    if storage_method == 's':
        pkl = pickle.dumps(noise10)
#       print(pkl)
    elif storage_method == 'p':
        pkl = pickle.dump(noise10, open("noise.pkl", "wb"))
    else:
        noise10.dump_hdf5("noise.h5")

    if storage_method == 's':
        clone = pickle.loads(pkl)
        #print()
        #print(clone.pet)
    elif storage_method == 'p':
        clone = pickle.load(open("noise.pkl", "rb"))
        #print()
        #print(clone.pet)
    else:
        noise10.load_hdf5("noise.h5")
        #print()
        #print(noise10.pet)


# test function. not a unit test.. yet
if __name__ == '__main__':
    import timeit

    orbit_range = [3300,53300]
#    orbit_range = [30000,31000]
    db_fname = "noise_.h5"

    # NOTE: entire db will get truncated!
    f = h5py.File(db_fname, "w")
    f.close()

    noise10 = Noise(1.0)
    noise10.extract(orbit_range)
    noise10.finalize()
    noise10.dump_hdf5(db_fname)

    pixnr=120
    print(noise10.noise[:,pixnr])
    print(np.min(noise10.noise[:,pixnr]), np.max(noise10.noise[:,pixnr]))
    print(noise10.noise[:,pixnr].shape)

    noise05 = Noise(0.5)
    noise05.extract(orbit_range)
    noise05.finalize()
    noise05.dump_hdf5(db_fname)

    noise05 = Noise(0.125)
    noise05.extract(orbit_range)
    noise05.finalize()
    noise05.dump_hdf5(db_fname)

#   print(timeit.timeit("test_load_store()", setup="from __main__ import test_load_store", number=10))
