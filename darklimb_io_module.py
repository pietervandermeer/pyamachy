# -*- coding: iso-8859-1 -*-
import numpy, h5py

# TODO: orbit as a list, or just keep it a scalar?
def sdmf_read_rts_darklimb_record(orbit, calib_db, 
    darklimbperiod=False, eclipseperiod=False, execmean=False, 
    execstddev=False, exectime=False, intersection=False, rtslevel=False, 
    jumpratios=False, peakratios=False, selected_orbits=False, pixelrange=[0,1024]):
 
    dict = {}
    dict['status'] = -1

    n_pixels = pixelrange[1]-pixelrange[0]+1
    p0 = pixelrange[0]
    p1 = pixelrange[1]

    fid = h5py.File(calib_db, 'r')

    #
    # find orbits.. 
    #

    orbitlist  = fid['orbitList'][:]
    meta_did = fid["metaTable"]

    #
    # TODO: find ranges in idx array (faster hyperslab processing)
    # 

    print orbitlist, orbit
    idx = numpy.where(orbitlist == orbit)
    print 'idx=', idx
    print 'idx[0].size=', idx[0].size 
    if idx[0].size is 0:
        print "sdmf_read_rts_darklimb_record: no orbits?"
        fid.close()
        return dict
    idx = (orbitlist == orbit)

    #
    # meta table
    #

    dict['mtbl'] = meta_did[idx]

    selected_orbits = orbitlist

    #
    # 2D arrays
    #

    if darklimbperiod:
        dict['darklimbperiod'] = fid['darkLimbPeriod'][:,idx]
    if eclipseperiod:
        dict['eclipseperiod']  = fid['eclipsePeriod'][:,idx]
    if intersection:
        dict['intersection']   = fid['intersection'][:,idx]
    if jumpratios:
        dict['jumpratios']     = fid['jumpRatios'][:,idx]
    if peakratios:
        dict['peakratios']     = fid['peakRatios'][:,idx]

    #
    # execution 3D arrays : executions x pixels x orbits  
    #

    if execmean:
        dict['execmean']   = fid["execMean"][:,:,idx]
    if execstddev:
        dict['execstddev'] = fid["execStddev"][:,:,idx]

    #
    # mean execution times : 2D : executions x orbits
    #

    if exectime:
        dict['exectime']   = fid["execTime"][:,idx]

    #
    # rts level : 3D : peaks x pixels x orbits
    # 

    if rtslevel:
        dict["rtslevel"]   = fid["rtsLevel"][:,:,idx]

    fid.close()
    dict['status'] = 0
    return dict
