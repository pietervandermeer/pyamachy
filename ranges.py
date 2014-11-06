from __future__ import print_function, division

def remove_overlap(ranges):
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

# test if range is fully contained in ranges
def is_in_range(ranges, range):
    for start, stop in sorted(ranges):
        if start <= range[0] and stop >= range[1]:
            return True
    return False

# unit test
if __name__ == "__main__":

    success = True;

    l = [(3,4),(0,1),(2,5)]
    new = remove_overlap(l)
    expected = [(0,1),(2,5)]
    if (new != expected):
        print("error!", new, " is not expected:", expected)
        success = False

    l.append((5,6))
    new = remove_overlap(l)
    expected = [(0,1),(2,6)]
    if (new != expected):
        print("error!", new, " is not expected:", expected)
        success = False

    l.append((1,2))
    new = remove_overlap(l)
    expected = [(0,6)]
    if (new != expected):
        print("error!", new, " is not expected:", expected)
        success = False

    l.append((10,15))
    new = remove_overlap(l)
    expected = [(0,6),(10,15)]
    if (new != expected):
        print("error!", new, " is not expected:", expected)
        success = False

    # get rid of all overlap and use this for range testing
    l = remove_overlap(l)

    # positive test
    test_range = (0,1)
    if not is_in_range(l, test_range):
        print("error!", test_range, " should be in:", l)
        success = False

    # negative test
    test_range = (11,16)
    if is_in_range(l, test_range):
        print("error!", test_range, " should not be in:", l)
        success = False

    print("unit test result =", success)
