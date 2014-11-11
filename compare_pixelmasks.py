from __future__ import print_function, division

import numpy as np

class Mask:

    def __init__(self):
        self.reset_pixel_window()
        return

    def load_ascii(self, orbit):
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

    # ASCII version of pixelquality 
    def load_ascii_quality(self, fname):
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
            return []

    def get_new_alive(self, new):
        tmp = self.mask.astype('i') - new.mask
        channel_indices = np.where(tmp == 1)[0]
        if channel_indices.size > 0:
            window = (channel_indices >= self.win_start) & (channel_indices <= self.win_end)
            return channel_indices[window]
        else:
            return []

    def apply_pixel_window(self, win_start, win_end):
        self.win_start = win_start
        self.win_end = win_end

    def reset_pixel_window(self):
        self.win_start = 0
        self.win_end = 1024

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

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    #print_new_dead_new_alive()
    print_dead_quality()

    # m2 = Mask()
    # m2.load_ascii(49245)
    # m3 = Mask()
    # m3.load_ascii(49259)

    # d1 = m1.diff(m2).astype('i')
    # d2 = m2.diff(m3).astype('i')
    # print(m1.mask)
    # print(m2.mask)
    # print(np.where(d1)[0])
    # print(np.where(d2)[0])

    # x = range(1024)
    # plt.cla()
    # plt.ticklabel_format(useOffset=False)
    # plt.ylim([-.1,1.1])
    # plt.plot(x, d1, ls='none', marker='o', label='diff1')
    # plt.plot(x, d2, ls='none', marker='o', label='diff2')
    # plt.legend(loc='best')
    # plt.show()
