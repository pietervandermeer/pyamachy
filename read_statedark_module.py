# -*- coding: iso-8859-1 -*-
#
# COPYRIGHT (c) 2011 SRON (pieter.van.der.meer@sron.nl)
#
#   This is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License, version 2, as
#   published by the Free Software Foundation.
#
#   The software is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, 
#   Boston, MA  02111-1307, USA.
#
import numpy
import h5py

#+
# NAME:
#       sdmf_read_statedark
#
# PURPOSE:
#       read rows of given orbit from statedark metaTable. In
#       addition, one can obtain other data parameters.
#
# CATEGORY:
#       SDMF - SCIA calibration
#
# CALLING SEQUENCE:
#       ToDo...
#
# INPUTS:
#       absorbit:    (absolute) Orbit number, long scalar or 2 element lonarr 
#                    (=orbit range)
#       state_id:    state id (integer in range [1..70])
#       calib_db:    path to the SDMF statedark database
#
# RETURNS:
#       dict:        dictionary of all requested data + metatable
#
# KEYWORD PARAMETERS:
#       status:      flag to return named variable with error flag (0 = ok)
#       smoothwidth: flag to return smoothing width (float, 1 per orbit)
#       stddev:      flag to return standard deviation
#       npeaks:      flag to return number of peaks in histogram
#       darkmedian:  flag to return median of dark signal
#       darksignal:  flag to return dark signal values at the two peaks
#       peakcounts:  flag to return histogram filling at the two peaks
#
# PROCEDURE:
#   ToDo...
#
# EXAMPLE:
#       data = sdmf_read_statedark([10000,10001], mtbl, 8, 
#           calib_db='/SCIA/share/SDMF/3.1/sdmf_statedark_ch6+.h5'),
#           stddev=True, npeaks=True)
#
# NOTES:
#       converted from IDL code.
#
# MODIFICATION HISTORY:
#   Written by: Pieter van der Meer (SRON), Oct 2011
#-
def sdmf_read_statedark(orbitrange, state_id, calib_db, 
                        stddev=False, npeaks=False, darkmedian=False,
                        darksignal=False, peakcounts=False, smoothwidth=False):

    if len(orbitrange) is not 2:
        print('sdmf_read_statedark_: orbitrange should have 2 elements')
        return

    dict = {}

    #
    # obtain indices to entries of the requested orbit and state
    #

    fid = h5py.File(calib_db, 'r') # generates an exception
    gid = fid["State_"+str('%02d'%state_id)] # generates an exception
    mt_did = gid['metaTable']
    orbitlist = (gid['orbitList'])[:]
    print('orbitrange=', orbitrange)
#    metaindx = numpy.where((orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1]))
    metaindx = (orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1])
    print(metaindx)
    if metaindx[0].size is 0:
        print('sdmf_read_statedark_: orbit range not present in database')
        dict['status'] = -1
        return dict

    mtbl = mt_did[metaindx]
    dict['mtbl'] = mtbl

    if darksignal:
        #darksignal = fltarr(1024,2,count)
        ds_did     = gid['DarkSignal']
        darksignal = ds_did[:,:,metaindx]
        dict['DarkSignal'] = darksignal

    if smoothwidth:
        #smoothwidth = fltarr(1024,count)
        sw_did     = gid['SmoothWidth']
        smoothwidth = sw_did[:,metaindx]
        dict['SmoothWidth'] = darksignal

    if darkmedian:
        #darkmedian = fltarr(1024,count)
        dm_did     = gid['DarkMedian']
        darkmedian = dm_did[:,metaindx]
        dict['DarkMedian'] = darksignal

    if stddev:
        #stddev     = fltarr(1024,2,count)
        sd_did     = gid['StdDev']
        stddev     = sd_did[:,:,metaindx]
        dict['StdDev'] = darksignal

    if npeaks:
        #npeaks     = bytarr(1024,count)
        np_did     = gid['Npeaks']
        dict['Npeaks'] = np_did[:,metaindx]

    if npeaks:
        #peakcounts = uintarr(1024,2,count)
        pc_did     = gid['PeakHeights']
        dict['PeakHeights'] = pc_did[:,:,metaindx]

    #
    # close file and return status and data
    #

    fid.close()
    dict['status'] = 0
    return(dict)
