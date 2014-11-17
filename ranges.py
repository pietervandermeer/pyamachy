"""
package that keeps track of ranges, can remove overlap between them and also merge adjacent ranges.
"""

from __future__ import print_function, division

import unittest

def remove_overlap(ranges):
    """
    Remove overlap, f.e. [(0,2),(2,3)] -> [(0,3)]
    """
    result = []
    current_start = -1
    current_stop = -1 

    for start, stop in sorted(ranges):
        if start > current_stop:
            # this segment starts after the last segment stops
            # just add a new segment
            result.append( (start, stop) )
            current_start, current_stop = start, stop
        else:
            # current_start already guaranteed to be lower
            current_stop = max(current_stop, stop)
            # segements overlap, replace
            result[-1] = (current_start, current_stop)

    return result

def merge_ranges(ranges):
    """
    Merge adjacent ranges, f.e. [(0,2),(3,5)] -> [(0,5)]
    """
    result = []
    current_start = -2 # needs to be 2 away from minimum
    current_stop = -2 # needs to be 2 away from minimum
    merged = ()

    for start, stop in sorted(ranges):
        if start == current_stop + 1:
            merged = (current_start, stop)
            current_start, current_stop = merged
        elif merged != ():
            result.append(merged)
            merged = ()
            current_start, current_stop = start, stop
        else:
            if current_start >= 0:
                result.append((current_start,current_stop))
            current_start, current_stop = start, stop

    if merged != ():
        result.append(merged)
    else:
        result.append((current_start, current_stop))
    return result

# test if range is fully contained in ranges
def is_in_range(ranges, range):
    for start, stop in sorted(ranges):
        print("is_in_range", range, start, stop)
        if start <= range[0] and stop >= range[1]:
            return True
    print(" not found")
    return False

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
