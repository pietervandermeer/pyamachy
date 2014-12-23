#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from ranges import remove_overlap, merge_ranges, is_in_range

class RangesTestCase(unittest.TestCase):

    def setUp(self):
        return

    def test_1(self):
        """ real-life test case """
        l = [(43360,43361)]
        l.append((43362,43370))
        l = remove_overlap(l)
        l = merge_ranges(l)
        expected = [(43360,43370)]
        self.assertEqual(l, expected)

        # l = [(43360,43361)]
        # l.append([43362,43370])
        # l = remove_overlap(l)
        # l = merge_ranges(l)
        # expected = [(43360,43370)]
        # self.assertEqual(l, expected)
        return

    def test_merge(self):
        l = [(0,2),(3,5)]
        new = merge_ranges(l)
        expected = [(0,5)]
        self.assertEqual(new, expected)

        l = [(0,2),(3,5),(7,10)]
        new = merge_ranges(l)
        expected = [(0,5),(7,10)]
        self.assertEqual(new, expected)

        l = [(0,2),(4,6),(7,10)]
        new = merge_ranges(l)
        expected = [(0,2),(4,10)]
        self.assertEqual(new, expected)

        l = [(0,2),(4,6),(7,10)]
        new = merge_ranges(l)
        expected = [(0,2),(4,10)]
        self.assertEqual(new, expected)

        l = [(0,2),(4,6),(7,10),(12,13)]
        new = merge_ranges(l)
        expected = [(0,2),(4,10),(12,13)]
        self.assertEqual(new, expected)

        l = [(0,2),(4,5)]
        new = merge_ranges(l)
        expected = [(0,2),(4,5)]
        self.assertEqual(new, expected)
        return

    def test_remove_overlap(self):
        l = [(3,4),(0,1),(2,5)]
        new = remove_overlap(l)
        expected = [(0,1),(2,5)]
        self.assertEqual(new, expected)

        l.append((5,6))
        new = remove_overlap(l)
        expected = [(0,1),(2,6)]
        self.assertEqual(new, expected)

        l.append((1,2))
        new = remove_overlap(l)
        expected = [(0,6)]
        self.assertEqual(new, expected)

        l.append((10,15))
        new = remove_overlap(l)
        expected = [(0,6),(10,15)]
        self.assertEqual(new, expected)

        # get rid of all overlap and use this for range testing
        l = remove_overlap(l)

        # positive test
        test_range = (0,1)
        self.assertTrue(is_in_range(l, test_range))

        # negative test
        test_range = (11,16)
        self.assertFalse(is_in_range(l, test_range))

        return

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
