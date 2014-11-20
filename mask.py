from __future__ import print_function, division

import numpy as np
import h5py

class Mask:

    def __init__(self):
        self.reset_pixel_window()
        return

    def load_ascii(self, orbit):
        # smooth mask
        fname = "/SCIA/SDMF30/Smoothmask/ASCII/"+str(orbit)+".mask"
        with open(fname) as f:
            content = f.readlines()

            #
            # skip header
            #

            content = content[6:]

            #
            # read in pixels
            #

            mask = np.zeros([8192], dtype=np.bool) 
            i = 0
            for line in content:
                mask[i] = bool(int(line.rstrip('\n')))
                i += 1

        #
        # store only channel 8
        #

        self.mask = mask[7*1024:8*1024]

        return

    def load_ascii_quality(self, fname):
        """
        load ASCII version of pixelquality 
        """

        with open(fname) as f:
            content = f.readlines()

            #
            # read in pixels
            #

            mask = np.zeros([1024], dtype=np.bool) 
            i = 0
            for line in content:
                mask[i] = bool(int(float(line.rstrip('\n'))))
                i += 1

        self.mask = mask

        return

    def load_sdmf30_crit(self, orbit, crit_name, smooth=False):
        """
        load sdmf 3.0 pixelmask's noise criterion  
        """
        fname = "/SCIA/SDMF30/sdmf_pixelmask.h5"
        fid = h5py.File(fname, "r")
        if smooth:
            grpname = "smoothMask/"
        else:
            grpname = "orbitalMask/"
        ds_orbits = fid[grpname+"orbitList"]

        orbits = ds_orbits[:]
        #print(orbits)

        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")

        i = np.argmax(idx)
        #print("i=",i)

        dset = fid[grpname+crit_name]
        self.mask = np.array(dset[7*1024:,i], dtype=np.bool).flatten()
        return

    def load_sdmf32_figure(self, orbit, fig_name):
        """
        load sdmf 3.2 pixelmask's figure (floats) and thresholds them to bools 
        """
        fname = "sdmf_pyxelmask.h5"
        fid = h5py.File(fname, "r")
        grpname = "" #orbitalMask/" # maybe in the future
        ds_orbits = fid[grpname+"orbits"]

        orbits = ds_orbits[:]
        #print(orbits)

        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")

        i = np.argmax(idx)
        print("i=",i)

        dset = fid[grpname+fig_name]
        fig = dset[i,:]
        # all nan's should get mask = 0 (not flagged), for good comparison with 3.0.. 
        idx = np.isfinite(fig)
        self.mask = np.zeros(1024, dtype=np.bool)
        self.mask[idx] = fig[idx] < 0.1 # quite arbitrary threshold.. 
        return

    def load_sdmf32_mask(self, orbit, crit_name):
        """
        load sdmf 3.2 pixelmask's mask (bools)
        """
        fname = "sdmf_pyxelmask.h5"
        fid = h5py.File(fname, "r")
        grpname = "" #orbitalMask/" # maybe in the future
        ds_orbits = fid[grpname+"orbits"]

        orbits = ds_orbits[:]
        #print(orbits)

        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")

        i = np.argmax(idx)
        print("i=",i)

        dset = fid[grpname+crit_name]
        self.mask = dset[i,:]
        return

    def diff(self, new):
        tmp = np.logical_xor(self.mask, new.mask)
        channel_indices = np.where(tmp)[0] 
        if channel_indices.size > 0:
            window = (channel_indices >= self.win_start) & (channel_indices <= self.win_end)
            return channel_indices[window]
        else:
            return []

    def get_new_dead(self, new):
        tmp = new.mask.astype('i') - self.mask
        channel_indices = np.where(tmp == 1)[0]
        if channel_indices.size > 0:
            window = (channel_indices >= self.win_start) & (channel_indices <= self.win_end)
            return channel_indices[window]
        else:
            return np.array([])

    def get_new_alive(self, new):
        tmp = self.mask.astype('i') - new.mask
        channel_indices = np.where(tmp == 1)[0]
        if channel_indices.size > 0:
            window = (channel_indices >= self.win_start) & (channel_indices <= self.win_end)
            return channel_indices[window]
        else:
            return np.array([])

    def apply_pixel_window(self, win_start, win_end):
        self.win_start = win_start
        self.win_end = win_end

    def reset_pixel_window(self):
        self.win_start = 0
        self.win_end = 1024
