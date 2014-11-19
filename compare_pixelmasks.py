from __future__ import print_function, division

import numpy as np
import h5py

from sciamachy_module import get_closest_state_exec, NoiseModel, petcorr, get_darkstateid

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

    def load_sdmf30_crit(self, orbit, crit_name):
        """
        load sdmf 3.0 pixelmask's noise criterion  
        """
        fname = "/SCIA/SDMF30/sdmf_pixelmask.h5"
        fid = h5py.File(fname, "r")
        grpname = "orbitalMask/"
        ds_orbits = fid[grpname+"orbitList"]

        orbits = ds_orbits[:]
        #print(orbits)

        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")

        i = np.argmax(idx)
        print("i=",i)

        dset = fid[grpname+crit_name]
        self.mask = np.array(dset[7*1024:,i], dtype=np.bool)
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

def compare_combined_30_vs_32(orbit):
    m30 = Mask()
    m30.load_sdmf30_crit(orbit, "combined")

    m32 = Mask()
    m32.load_sdmf32_mask(orbit, "combinedFlag")

    # for i in range(1024):
    #     print(i, m30.mask[i], m32.mask[i])

    newdead = m30.get_new_dead(m32)
    newalive = m30.get_new_alive(m32)
    print("3.0->3.2 dead:", newdead.size, newdead)
    print("3.0->3.2 alive:", newalive.size, newalive)
    return

def compare_noise_30_vs_32(orbit):
    m30 = Mask()
    m30_combined = Mask()
    m30_combined.load_ascii(orbit)
    m30.load_sdmf30_crit(orbit, "noise")

    m32 = Mask()
    m32.load_sdmf32_figure(orbit, "noise")

    newdead = m30.get_new_dead(m32)
    newalive = m30.get_new_alive(m32)
    print("new noisy:", newdead.size, newdead)
    print("new noisy but already bad", np.sum(np.in1d(newdead, np.where(m30_combined.mask))))
    print("real new noisy", newdead[np.in1d(newdead, np.where(m30_combined.mask), invert=True)])

    print("new not-noisy:", newalive.size, newalive)
    print("new not-noisy but still bad", np.sum(np.in1d(newalive, np.where(m30_combined.mask))))
    print("real new not-noisy:", newalive[np.in1d(newalive, np.where(m30_combined.mask), invert=True)])
    return

def compare_res_30_vs_32(orbit):
    m30 = Mask()
    m30_combined = Mask()
    m30_combined.load_ascii(orbit)
    m30.load_sdmf30_crit(orbit, "residual")

    m32 = Mask()
    m32.load_sdmf32_figure(orbit, "darkResidual")

    newdead = m30.get_new_dead(m32)
    newalive = m30.get_new_alive(m32)
    print("new resy:", newdead.size, newdead)
    print("new resy but already bad", np.sum(np.in1d(newdead, np.where(m30_combined.mask))))
    print("real new resy", newdead[np.in1d(newdead, np.where(m30_combined.mask), invert=True)])

    print("new not-resy:", newalive.size, newalive)
    print("new not-resy but still bad", np.sum(np.in1d(newalive, np.where(m30_combined.mask))))
    print("real new not-resy:", newalive[np.in1d(newalive, np.where(m30_combined.mask), invert=True)])
    return

def compare_error_30_vs_32(orbit):
    m30 = Mask()
    m30_combined = Mask()
    m30_combined.load_ascii(orbit)
    m30.load_sdmf30_crit(orbit, "residual")

    m32 = Mask()
    m32.load_sdmf32_figure(orbit, "darkError")

    newdead = m30.get_new_dead(m32)
    newalive = m30.get_new_alive(m32)
    print("new errory:", newdead.size, newdead)
    print("new errory but already bad", np.sum(np.in1d(newdead, np.where(m30_combined.mask))))
    print("real new errory", newdead[np.in1d(newdead, np.where(m30_combined.mask), invert=True)])

    print("new not-errory:", newalive.size, newalive)
    print("new not-errory but still bad", np.sum(np.in1d(newalive, np.where(m30_combined.mask))))
    print("real new not-errory:", newalive[np.in1d(newalive, np.where(m30_combined.mask), invert=True)])
    return

def print_dead_quality():
    orbit = 50000
    m = Mask()
    m.load_ascii_quality("qualities/"+str(orbit)+".txt")
    print(np.where(m.mask))
    print(np.sum(m.mask))

    wlsmask = np.zeros(1024, dtype=np.bool)
    wls_idx = np.array([143, 196, 219, 220, 222, 254, 338, 398, 399, 405, 406, 494, 502, 577, 601, 609, 624, 635, 678, 707, 779, 802, 868, 881, 902])
    sun_idx = np.array([218, 219, 220, 609, 779, 883, 902, 1008])
    wlsmask[wls_idx] = True
    new_alive = m.mask.astype('i') - wlsmask
    print("new alive = ", np.where(new_alive == 1), np.where(new_alive == 1)[0].size)
    new_dead = wlsmask.astype('i') - m.mask
    print("new dead = ", np.where(new_dead  == 1), np.where(new_dead == 1)[0].size)
    print(m.mask)
    print(new_alive)
    print(new_dead)
    return

def print_new_dead_new_alive():
    start_orbit = 21847 
    end_orbit = 21890

#    start_orbit = 22820 
#    end_orbit = 22890

#    start_orbit = 24023
#    end_orbit = 24051

#    start_orbit = 27229 
#    end_orbit = 27258

#    start_orbit = 32726 
#    end_orbit = 32870

#    start_orbit = 49230
#    end_orbit = 49259

    m1 = Mask()
    m1.apply_pixel_window(394, 620)
    m1.load_ascii(start_orbit)

    m2 = Mask()
    for orbit in range(start_orbit, end_orbit):
        m2.load_ascii(orbit)
        print(orbit, "new dead:", m1.get_new_dead(m2), "new alive:", m1.get_new_alive(m2))

    return

def print_criteria32(orbit, pixnr):
    """ 
    Print flagging criteria for SDMF3.2.
    """
    pixnr -= 7*1024
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

    ds_combinedf = fid[grpname+"combinedFlag"]
    print("flag",ds_combinedf[i,pixnr])
    ds_combined = fid[grpname+"combined"]
    print("q",ds_combined[i,pixnr])
    ds_inv = fid[grpname+"invalid"]
    print("inv",ds_inv[i,pixnr])
    ds_dc = fid[grpname+"saturation"]
    print("sat",ds_dc[i,pixnr])
    ds_noise = fid[grpname+"noise"]
    print("noise",ds_noise[i,pixnr])
    ds_res = fid[grpname+"darkResidual"]
    print("res",ds_res[i,pixnr])
    ds_res = fid[grpname+"darkError"]
    print("error",ds_res[i,pixnr])
    ds_sun = fid[grpname+"sunResponse"]
    print("sun",ds_sun[i,pixnr])
    ds_wls = fid[grpname+"wlsResponse"]
    print("wls",ds_wls[i,pixnr])

    return    

def print_criteria30(orbit, pixnr):
    """ 
    Print flagging criteria for SDMF3.0.
    """
    m = Mask()

    #
    # first just print combined smoothmask
    #

    m.load_ascii(orbit)
    print(np.where(m.mask))

    fname = "/SCIA/SDMF30/sdmf_pixelmask.h5"
    fid = h5py.File(fname, "r")
    grpname = "orbitalMask/"
    ds_orbits = fid[grpname+"orbitList"]

    orbits = ds_orbits[:]
    #print(orbits)

    idx = orbit == orbits
    if np.sum(idx) == 0:
        raise Exception("orbit not in orbital mask")

    i = np.argmax(idx)
    print("i=",i)

    ds_combined = fid[grpname+"combined"]
    print("combined flag",ds_combined[pixnr,i])
    ds_ao = fid[grpname+"analogOffset"]
    print("ao",ds_ao[pixnr,i])
    ds_aoerr = fid[grpname+"analogOffsetError"]
    print("aoerr",ds_aoerr[pixnr,i])
    ds_dc = fid[grpname+"darkCurrent"]
    print("dc",ds_dc[pixnr,i])
    ds_dcerr = fid[grpname+"darkCurrentError"]
    print("dcerr",ds_dcerr[pixnr,i])
    ds_chi = fid[grpname+"chiSquare"]
    print("chi",ds_chi[pixnr,i])
    ds_noise = fid[grpname+"noise"]
    print("noise",ds_noise[pixnr,i])
    ds_inv = fid[grpname+"invalid"]
    print("inv",ds_inv[pixnr,i])
    ds_ppg = fid[grpname+"pixelGain"]
    print("ppg",ds_ppg[pixnr,i])
    ds_res = fid[grpname+"residual"]
    print("res",ds_res[pixnr,i])
    ds_sun = fid[grpname+"sunResponse"]
    print("sun",ds_sun[pixnr,i])
    ds_wls = fid[grpname+"wlsResponse"]
    print("wls",ds_wls[pixnr,i])

    return

def print_nr_ppg(orbit):

    dictwls = get_closest_state_exec(orbit, 61, "/SCIA/SDMF31/sdmf_extract_calib.h5", readoutMean=True)
    orbit_range = dictwls["orbit_range"]
    orbit = orbit_range[0]
    print("closest wls orbit=", orbit)

    m = Mask()

    #
    # first just print combined smoothmask
    #

    m.load_ascii(orbit)
    print(np.where(m.mask))

    fname = "/SCIA/SDMF30/sdmf_pixelmask.h5"
    fid = h5py.File(fname, "r")
    grpname = "orbitalMask/"
    ds_orbits = fid[grpname+"orbitList"]

    orbits = ds_orbits[:]
    #print(orbits)

    idx = orbit == orbits
    if np.sum(idx) == 0:
        raise Exception("orbit not in SDMF3.2 orbital mask")

    i = np.argmax(idx)
    print("i=",i)

    ds_ppg = fid[grpname+"pixelGain"]
    ppg_chan8 = ds_ppg[7*1024:,i]
    print("ppg", np.sum(ppg_chan8), np.where(ppg_chan8))

    return

def print_dark(orbit, pixnr):

    #import matplotlib.pyplot as plt

    fname = "/SCIA/SDMF30/sdmf_dark.h5"
    fid = h5py.File(fname, "r")
    ds_orbits = fid["orbitList"]

    orbits = ds_orbits[:]
    #print(orbits)

    idx = orbit == orbits
    if np.sum(idx) == 0:
        raise Exception("orbit not in SDMF3.0 orbital mask")

    i = np.argmax(idx)
    print("i=",i)
    ds_ao = fid["analogOffset"]
    ds_aoerr = fid["analogOffsetError"]
    ds_dc = fid["darkCurrent"]
    ds_dcerr = fid["darkCurrentError"]
    ds_chi = fid["chiSquareFit"]

    print("dark=",ds_dc[pixnr,i])
    print("ao=",ds_ao[pixnr,i])
    print("chi=",ds_chi[pixnr,i])

    # plt.cla()
    # plt.plot(ds_dc[pixnr,i])
    # plt.show()

    return

def print_noise30(orbit, pixnr):
    pets = 0.125, 0.5, 1.0
    for pet in pets:
        stateid = get_darkstateid(pet, orbit)
        fname = "/SCIA/SDMF30/sdmf_extract_calib.h5"
        fid = h5py.File(fname, "r")

        grpname = "State_"+str('%02d'%stateid)+"/"
        ds_orbits = fid[grpname+"orbitList"]
        orbits = ds_orbits[:]
        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")
        i = np.argmax(idx)
        print("i=",i)
        ds_noise = fid[grpname+"readoutNoise"]
        print("SDMF3.0 noise("+str(pet)+")=", ds_noise[pixnr,i])

    return

def print_noise(orbit, pixnr):
    fname = "/SCIA/SDMF31/pieter/noise.h5"
    fid = h5py.File(fname, "r")
    pets = 0.125, 0.5, 1.0

    for pet in pets:
        ds_orbits = fid["pet"+str(pet)+"/orbits"]
        orbits = ds_orbits[:]
        idx = orbit == orbits
        if np.sum(idx) == 0:
            raise Exception("orbit not in orbital mask")
        i = np.argmax(idx)
        print("i=",i)
        ds_noise = fid["pet"+str(pet)+"/noise"]
        print("SDMF3.2 noise("+str(pet)+")=", ds_noise[i, pixnr-7*1024])

    return

def print_vardark(orbit, pixnr):
    fname = "/SCIA/SDMF31/pieter/vardark_long.h5"
    fid = h5py.File(fname, "r")

    ds_orbits = fid["dim_orbit"]
    orbits = ds_orbits[:]
    #print(orbits)
    idx = orbit == orbits
    if np.sum(idx) == 0:
        raise Exception("orbit not in orbital mask")
    i = np.argmax(idx)
    print("i=",i)
    ds_dark = fid["varDark"]
    print("dark (sdmf3.2)=", ds_dark[i, 0, pixnr-7*1024])
    ds_un = fid["uncertainties"]
    print("uncertainty=", ds_un[i, pixnr-7*1024])

    return

def print_30_vs_32_fit_quality(orbit):
    """ 
    Print flagging criteria for SDMF3.0 vs SDMF3.2 for each pixel.
    """

    fname = "/SCIA/SDMF30/sdmf_pixelmask.h5"
    fid30 = h5py.File(fname, "r")
    grpname = "orbitalMask/"
    ds_orbits = fid30[grpname+"orbitList"]
    orbits = ds_orbits[:]
    idx = orbit == orbits
    if np.sum(idx) == 0:
        raise Exception("orbit not in orbital mask")
    i30 = np.argmax(idx)

    ds_combined30 = fid30[grpname+"combined"]
    ds_chi30 = fid30[grpname+"chiSquare"]
    ds_noise30 = fid30[grpname+"noise"]
    ds_inv30 = fid30[grpname+"invalid"]
    ds_res30 = fid30[grpname+"residual"]

    fname = "sdmf_pyxelmask.h5"
    fid32 = h5py.File(fname, "r")
    grpname = ""
    ds_orbits = fid32[grpname+"orbits"]
    orbits = ds_orbits[:]
    idx = orbit == orbits
    if np.sum(idx) == 0:
        raise Exception("orbit not in orbital mask")
    i32 = np.argmax(idx)

    ds_combined32 = fid32[grpname+"combinedFlag"]
    ds_res32 = fid32[grpname+"darkResidual"]
    ds_error32 = fid32[grpname+"darkError"]
    ds_noise32 = fid32[grpname+"noise"]

    print("pix\tcom\t\tchi\t\tnoi\t\tinv\t\tres\t\tc32\t\tr32\t\te32\t\tn32")
    for pixnr in range(1024):
        print(pixnr, "\t\t", 
              ds_combined30[pixnr+7*1024,i30], "\t\t", 
              ds_chi30[pixnr+7*1024,i30], "\t\t", 
              ds_noise30[pixnr+7*1024,i30], "\t\t", 
              ds_inv30[pixnr+7*1024,i30], "\t\t", 
              ds_res30[pixnr+7*1024,i30], "\t\t",
              ds_combined32[i32,pixnr], "\t\t", 
              ds_res32[i32,pixnr], "\t\t", 
              ds_error32[i32,pixnr], "\t\t", 
              ds_noise32[i32,pixnr])

    return

def print_noisemodel(pixnr):
    nm = NoiseModel()
    print("noise model @ 1.0s:", nm.compute(pixnr, 1.0-petcorr, 6000.))
    print("noise model @ 0.5s:", nm.compute(pixnr, 0.5-petcorr, 6000.))
    print("noise model @ 0.125s:", nm.compute(pixnr, 0.125-petcorr, 6000.))
    return

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    #print_new_dead_new_alive()
    #print_dead_quality()
    orbit = 43000
    pixnr = 450+7*1024

    # print_nr_ppg(orbit)
    #print_criteria30(orbit, pixnr)
    #print_dark(orbit, pixnr)
    #print_noise(orbit, pixnr)
    print_vardark(orbit, pixnr)
    print_criteria32(orbit, pixnr)
    # compare_noise_30_vs_32(orbit)
    # compare_res_30_vs_32(orbit)
    # compare_error_30_vs_32(orbit)
    print_noise30(orbit, pixnr)
    compare_combined_30_vs_32(orbit)
    #print_30_vs_32_fit_quality(orbit)

    print_noisemodel(pixnr)
