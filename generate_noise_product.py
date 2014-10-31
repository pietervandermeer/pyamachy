from __future__ import print_function, division

import cPickle as pickle # python 2 speed-up
import pickle
import numpy as np 
import h5py

from sciamachy_module import get_darkstateid, read_extracted_states_
from envisat import PhaseConverter

#-- globals --------------------------------------------------------------------

n_pix = 1024

# storage method. interesting.. pickle is incompatible with ctypes.. so never mind.. 
#storage_method = 's' # pickle string
#storage_method = 'p' # pickle binary file
storage_method = 'h' # hdf5

#-- functions ------------------------------------------------------------------

class Noise:
    def __init__(self, pet):
        self.pet = pet
        self.calib_db = '/SCIA/SDMF31/sdmf_extract_calib.h5'
        self.pc = PhaseConverter()

    # extract noise data from sdmf_extract_calib database
    def extract(self, orbit_range):
        if orbit_range[0] < 43362 and orbit_range[1] >= 43362:
            pre_ocr43_sid = get_darkstateid(self.pet, 43361)
            post_ocr43_sid = get_darkstateid(self.pet, 43362)
            pre_range = [orbit_range[0], 43361]
            post_range = [43362, orbit_range[1]]
            d1 = read_extracted_states_(pre_range, pre_ocr43_sid, self.calib_db, readoutNoise=True, errorInTheMean=False, readoutCount=True)
            d2 = read_extracted_states_(post_range, post_ocr43_sid, self.calib_db, readoutNoise=True, errorInTheMean=False, readoutCount=True)
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
            d = read_extracted_states_(orbit_range, stateid, self.calib_db, readoutNoise=True, errorInTheMean=False, readoutCount=True)
            self.count = d['readoutCount'][:]
            self.mjds = d['mtbl'][:]['julianDay']
            self.noise = d['readoutNoise'][:,:]
        ephases, self.orbits = self.pc.get_phase(self.mjds, getOrbits=True)

    # averages noise figures so that they are per orbit
    def finalize(self):
        # sort.. to group state executions by orbit
        idx = np.argsort(self.orbits)
        self.orbits = self.orbits[idx]
        self.noise = self.noise[idx,:]

        #
        # average the noise figure for each orbit
        #

        # get position and group size
        firsts = np.nonzero(np.diff(self.orbits))[0] + 1
        firsts = np.insert(firsts, 0, 0)
        diffs = np.diff(firsts)
        counts = np.insert(diffs, -1, self.orbits.size - firsts[-1])

        # create output array (where all groups are collapsed into single columns)
        n_orbits = counts.size
        noise_out = np.empty([n_orbits, n_pix])

        # average each orbit
        for i in range(n_orbits):
            count = counts[i]
            first = firsts[i]
            sdev = self.noise[first:first+count, :]  # stddev
            count = np.sum(self.count[first:first+count, :], axis=0)  # total count of all measurements
            noise_out[i,:] = np.sqrt(np.mean(np.square(sdev), axis=0) / count) # root(mean_square / count) -> error in the mean

        self.noise = noise_out

    def get_groupname(self):
        return "pet"+str(self.pet)

    def dump_hdf5(self, fname):
        f = h5py.File(fname, "a")
        grpn = self.get_groupname()
        n_dset = f.create_dataset(grpn+"/noise", self.noise.shape, dtype='f')
        n_dset[:,:] = self.noise
        o_dset = f.create_dataset(grpn+"/orbits", self.orbits.shape, dtype='i')
        o_dset[:] = self.orbits
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
    db_fname = "noise.h5"

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
