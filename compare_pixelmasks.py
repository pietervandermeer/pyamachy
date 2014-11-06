from __future__ import print_function, division

import numpy as np

class Mask:

    def __init__(self):
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
                #print(i, bool(int(line.rstrip('\n'))), int(line.rstrip('\n')), line.rstrip('\n'))
                mask[i] = bool(int(line.rstrip('\n')))
                i += 1

        #
        # store only channel 8
        #

        self.mask = mask[7*1024:8*1024]
        return

    def diff(self, m2):
        return np.logical_xor(self.mask, m2.mask) 

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

    start_orbit = 27229 
    end_orbit = 27258

    m1 = Mask()
    m1.load_ascii(start_orbit)

    m2 = Mask()
    for orbit in range(start_orbit, end_orbit):
        m2.load_ascii(orbit)
        diff = m1.diff(m2)
        print(orbit, np.where(diff)[0])

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
