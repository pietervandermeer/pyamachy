#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from vardark_module import AllDarks

pixnr= 535

ad = AllDarks([0.125, 0.5, 1.0])

orbit_range = [43300,43400]
n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = ad.get_range(orbit_range)

plt.cla()
idx_0125 = (np.where(pet <=   .125))[0]
idx_0500 = (np.where((pet <=  .500) & (pet > .125)))[0]
idx_1000 = (np.where((pet <= 1.000) & (pet > .500)))[0]
plt.plot(ephases[idx_0125], all_readouts[idx_0125,pixnr], 'bo')
plt.plot(ephases[idx_0500], all_readouts[idx_0500,pixnr], 'go')
plt.plot(ephases[idx_1000], all_readouts[idx_1000,pixnr], 'ro')
print("0.125:", all_readouts[idx_0125,pixnr])
print("0.5:", all_readouts[idx_0500,pixnr])
print("1.0:", all_readouts[idx_1000,pixnr])
plt.show()

