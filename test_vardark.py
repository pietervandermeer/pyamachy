#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import unittest
import numpy as np
import numpy.testing as nptst

from sciamachy import petcorr
from vardark import AllDarks

class AllDarksTestCase(unittest.TestCase):

    def setUp(self):
        return

    def test_data_range(self):
        """ test if returned data is in correct orbit range """
        in_pets = [0.125, 0.5, 1.0]
        orbit_range = [15000, 15001]
        ad = AllDarks(in_pets)
        n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = ad.get_range(orbit_range)
        self.assertGreater(n_exec, 0)
        self.assertTrue(np.all(ephases >= orbit_range[0]))
        self.assertTrue(np.all(ephases <= orbit_range[1]))
        return

    def test_petcorr(self):
        """ test if pet timing correction is applied  """
        in_pets = [0.125, 0.5, 1.0]
        ad = AllDarks(in_pets)
        n_exec, all_state_phases, pet, coadd, all_readouts, all_sigmas, ephases = ad.get_range([15000,15001])
        self.assertGreater(n_exec, 0)
        pet_unique_uncorr = np.unique(np.sort(pet)) + petcorr
        nptst.assert_allclose(pet_unique_uncorr, in_pets, rtol=1e-5) # this tolerance should be ok for float32
        return

def visual_test():
    import matplotlib.pyplot as plt

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
    return

#-- main -------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
