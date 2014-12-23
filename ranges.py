#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
package that keeps track of ranges, can remove overlap between them and also merge adjacent ranges.
"""

from __future__ import print_function, division

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
        #print("is_in_range", range, start, stop)
        if start <= range[0] and stop >= range[1]:
            return True
    #print(" not found")
    return False
