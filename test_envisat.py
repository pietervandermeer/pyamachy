#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy as np
import numpy.testing as nptst

from envisat import PhaseConverter

#-- unit tests -----------------------------------------------------------------

class PhaseConverterTestCase(unittest.TestCase):

    def setUp(self):
        self.pc = PhaseConverter()
        return

    def test_fat1(self):
        """ 
        Just a randomly chosen factor acceptance test 
        """
        phases, orbits = self.pc.get_phase(np.array([1000.0,1000.1,1000.2,1000.3]), getOrbits=True)
        expected_phases = np.array([0.10914163, -0.45960116, -0.02825826, -0.59700102], dtype=np.float32)
        expected_orbits = [3004, 3006, 3007, 3009]
        self.assertEqual(orbits.tolist(), expected_orbits)
        nptst.assert_allclose(phases, expected_phases, rtol=1e-5) # this tolerance should be ok for float32
        return

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
