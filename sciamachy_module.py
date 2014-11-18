# -*- coding: iso-8859-1 -*-
#
# COPYRIGHT (c) 2014 SRON (pieter.van.der.meer@sron.nl)
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

from __future__ import print_function, division

import numpy as np
import h5py
import ConfigParser

"""
functions and classes related to SCIAMACHY operation and instrument
2011, 2014 - Pieter van der Meer - SRON - Netherlands Institute for Space Research
"""

channels = ['1','2','3','4','5','6','6+','7','8']
pixranges = [np.arange(1024), \
             1024+np.arange(1024), \
             1024*2+np.arange(1024), \
             1024*3+np.arange(1024), \
             1024*4+np.arange(1024), \
             1024*5+np.arange(795), \
             1024*5+795+np.arange(229), \
             1024*6+np.arange(1024), \
             1024*7+np.arange(1024), \
            ]
# exposure time error
petcorr = 1.18125e-3
# nr of pixels per channel
n_chanpix = 1024

# start orbits of the decontamination periods over the entire mission
decon_start_orbits = [4380, 5718, 6384, 7574, 9407, 12031, 14675, 35574]
# end orbits of the decontamination periods over the entire mission
decon_end_orbits = [4428, 5766, 6449, 7827, 9673, 12208, 14912, 35848]
# intervals between and after decontaminations
decon_intervals = [[4428,5718],[5766,6384],[6449,7574],[7827,9407],[9673,12031],[12208,14675],[14912,35574],[35848,52000]]

# TODO: this shouldn't be here. not directly related to the physical properties of the instrument
mask_criteria = ['combined','RTS','darkCurrentError','darkCurrentSat','invalid','residual']

class MemCorrector:
    
    def __init__(self):
        
        #
        # Parse config file, exit if unsuccessful
        #
        
        config_fname = open('default3.1.cfg')
        try:
            self.cfg = self.get_config(config_fname)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+config_fname+"'!")
            raise
        
        #
        # load memory correction table
        #
        
        fid = h5py.File(self.cfg['memcorr_fname'], 'r')
        dset = fid['MemTable']
        self.memtbl = dset[:]
        #print "memtbl.shape=", self.memtbl.shape
        fid.close()

    def correct(self, spectrum):
        for i in range(5):
            spec = (spectrum[i*1024:(i+1)*1024]).round().astype(int)
            spec = np.clip(spec,0,65535)
            spectrum[i*1024:(i+1)*1024] -= self.memtbl[i,spec]
        return spectrum

    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['memcorr_fname'] = parser.get('Global','memcorr_file')
        return dict

class NonlinCorrector:
    
    def __init__(self):
        
        #
        # Parse config file, exit if unsuccessful
        #
        
        config_fname = open('default3.1.cfg')
        try:
            self.cfg = self.get_config(config_fname)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+config_fname+"'!")
            raise
        
        #
        # load non-linearity correction tables
        #
        
        fid = h5py.File(self.cfg['nlcorr_fname'], 'r')
        dset = fid['CurveIndex']
        self.curveIndex = dset[:]
#        indx = curveIndex[nchan-1,ipix % 1024]
#        self.nlintbl = dset[indx,:]
        dset = fid['nLinTable']
        self.nlintbl = dset[:]
        #print "nlintable.shape", self.nlintbl.shape
        #print "curveIndex.shape", self.curveIndex.shape
        fid.close()

    def correct(self, spectrum):
        for i_chan in range(5,8):
            tabidx = self.curveIndex[i_chan,:]
            #print "tabidx.shape=", tabidx.shape
            #print tabidx
            spec = (spectrum[i_chan*1024:(i_chan+1)*1024]).round().astype(int)
            spec = np.clip(spec,0,65535)
            #print spec,spec.shape
            spectrum[i_chan*1024:(i_chan+1)*1024] -= self.nlintbl[tabidx,spec]
        return spectrum

    def correct_ch8(self, spectrum):
        tabidx = self.curveIndex[7,:]
        #print "tabidx.shape=", tabidx.shape
        #print tabidx
        spec = spectrum.round().astype(int)
        spec = np.clip(spec,0,65535)
        #print spec,spec.shape
        spectrum -= self.nlintbl[tabidx,spec]
        return spectrum

    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['nlcorr_fname'] = parser.get('Global','nlcorr_file')
        return dict

class pixelmask:
    def __init__(self):
        # SDMF 3.0 pixelmask: very useful indeed
        fmask = h5py.File('/SCIA/share/SDMF/3.0/sdmf_pixelmask.h5', 'r')
        gmask = fmask['smoothMask']
        self.orbits = (gmask['orbitList'])[:]
        self.msk = np.transpose(gmask['combined'][:,:])
        fmask.close()

class orbitfilter:
    """
    class which filters out special orbits. typically decontaminations and 
    monthly calibration orbits and such.
    """
    
    def __init__(self):
        
        #
        # Parse config file, exit if unsuccessful
        #
        
        config_fname = open('default3.1.cfg')
        try:
            self.cfg = self.get_config(config_fname)
        except ConfigParser.NoOptionError, ex:
            msg = "There was a missing option in the configuration file '"
            logging.exception(msg+config_fname+"'!")
            raise

        #
        # read orbit filter files into numpy arrays
        #

        # easy peasy
        full_monthlies_fname = self.cfg['db_dir']+self.cfg['monthlies_fname']
        # this one is tricky since it contains ranges along with HK info
        full_quality_fname   = self.cfg['db_dir']+self.cfg['quality_fname']
        self.monthlies = np.loadtxt(full_monthlies_fname, dtype=int)
        qua_file = open(full_quality_fname,'r')
        qua = qua_file.read()
        list1 = qua.split('\n')
        n_rows = len(list1)-1 # assuming there's a trailing \n
        qua_ranges = np.zeros((n_rows,2),dtype=int)
        #print qua_ranges.shape
        i=0
        for line in list1:
            columns = line.split()
            if len(columns) < 2:
                break
            #print columns[0:2]
            qua_ranges[i,0] = int(columns[0])
            qua_ranges[i,1] = int(columns[1])
            i+=1
        self.qualities = qua_ranges
        
    # load configuration from file
    def get_config(self, config_file):
        parser=ConfigParser.SafeConfigParser()
        parser.readfp(config_file)
        dict = {}
        dict['db_dir'] = parser.get('Global','masterdirectory')
        dict['quality_fname'] = parser.get('Global','quality_file')
        dict['monthlies_fname'] = parser.get('Global','monthlies_file')
        return dict

    # stamp out all the low-quality orbit ranges
    def get_quality_orbit_filter(self, orbits):
        mask = np.ones(orbits.size, dtype=bool)
        n_ranges = self.qualities.shape[0]
        for i in range(n_ranges):
            start = self.qualities[i,0]
            end = self.qualities[i,1]
            #loc_mask = (orbits < start) | (orbits > end)
            #mask &= loc_mask
            loc_idx = np.where((orbits >= start) & (orbits < end))
            mask[loc_idx] = False
        return mask

    # stamp out all the orbits around the monthly calibrations 
    # (without limb or nadir measurements)
    def get_monthly_orbit_filter(self, orbits):
        monthly_kernels = np.concatenate((self.monthlies,self.monthlies+2,
                                          self.monthlies+4,self.monthlies-2,
                                          self.monthlies-4))
        mask = np.in1d(orbits, monthly_kernels)
        mask = np.invert(mask)
        return mask

    # return closest monthly calibration orbit
    def get_closest_monthly(self, orbit):
        delta = np.abs(self.monthlies - orbit)
        idx = np.argmin(delta)
        if np.isscalar(idx):
            return self.monthlies[idx]
        else:
            return self.monthlies[idx[0]]

    # return closest monthly calibration orbit
    def get_previous_monthly(self, orbit):
        delta = self.monthlies - (orbit-1) # 2 orbits after eachother, so skip one.
        idx = np.where(delta < 0)
        if np.isscalar(idx):
            return self.monthlies[idx]
        elif idx[0].size == 0:
            # no earlier monthly found, return given orbit
            return orbit
        else:
            # return last monthly before given orbit
            return self.monthlies[idx[0][-1]]

    # return closest monthly calibration orbit
    def get_next_monthly(self, orbit):
        delta = self.monthlies - (orbit+1) # 2 orbits after eachother, so skip one. 
        idx = np.where(delta > 0)
        if np.isscalar(idx):
            #print("scal")
            return self.monthlies[idx]
        elif idx[0].size == 0:
            #print("00000")
            # no later monthly found, return given orbit
            return orbit
        else:
            #print("tuple:",self.monthlies[idx[0]])
            # return first monthly after given orbit
            return self.monthlies[idx[0][0]]

class OrbitRangeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# for given pet and absolute orbit number, return the the right dark state id
def get_darkstateid(pet, orbit):
    orbit = int(orbit)
    if orbit >= 1572 and orbit < 43362:
        if pet == 2:
            return 67
        elif pet == 1:
            return 8
        elif pet == 0.5:
            return 63
        elif pet == 0.125:
            return 26
        elif pet == 0.0625:
            return 46
        else:
            raise PetNotFoundError("pet not present between orbits 1572 and 43362")
    elif orbit >= 43362 and orbit < 60000:
        if pet == 1:
            return 67
        elif pet == 0.5:
            return 26
        elif pet == 0.375:
            return 63
        elif pet == 0.125:
            return 8
        elif pet == 0.0625:
            return 46
        else:
            raise PetNotFoundError("pet not present after orbit 43362")
    else:
        # very early orbit?
        raise OrbitRangeError("orbit "+str(orbit)+" out of range")

# reads extracted state executions from database
def read_extracted_states(orbitrange, state_id, calib_db, in_orbitlist=None, readoutMean=False, readoutNoise=False, orbitList=False):

    if in_orbitlist is None:
        if len(orbitrange) is not 2:
            raise Exception('orbitrange should have 2 elements')

    dict = {}

    #
    # obtain indices to entries of the requested orbit and state
    #

    fid = h5py.File(calib_db, 'r') # generates an exception
    gid = fid["State_"+str('%02d'%state_id)] # generates an exception
    mt_did = gid['metaTable']
    orbitlist = (gid['orbitList'])[:]
    #print(orbitlist)
    if in_orbitlist is not None:
        #metaindx = orbitlist = in_orbitlist
        metaindx = np.in1d(orbitlist, in_orbitlist)
    else:
        #print('orbitrange=', orbitrange)
        metaindx = (orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1])
        #print(metaindx, metaindx[0].size, np.sum(metaindx))

    if np.sum(metaindx) == 0:
        raise Exception('orbit range not present in database')

    mtbl = mt_did[metaindx]
    dict['mtbl'] = mtbl

    if readoutMean:
        ds_did      = gid['readoutMean']
        dict['readoutMean'] = ds_did[metaindx,:]

    # export error-in-the-mean instead of plain stddev.
    if readoutNoise:
        ds_did       = gid['readoutNoise']
        readoutNoise = ds_did[metaindx,:]
        ds_did       = gid['readoutCount']
        readoutNoise /= np.sqrt(ds_did[metaindx,:])
        readoutNoise *= 1.4826 # sigma = MAD * K 
        dict['readoutNoise'] = readoutNoise

    if orbitList:
        ds_did = gid['orbitList']
        dict['orbitList'] = ds_did[metaindx]

    #
    # find the pet and coadd (TODO: do for multiple orbits)
    #

    orbit = orbitrange[0]
    clusoff1 = [0,10,1014,1024]
    ds_did = gid['clusConf']
    clusconf = ds_did[:]
    orbit_start = clusconf['orbit'][:]
    cluspets = clusconf['pet'][:]
    cluscoad = clusconf['coaddf'][:]
    n_defchange = orbit_start.size 
    pet = np.empty(1024)
    coadd = np.empty(1024)
    if n_defchange is 1:
        orbit_end = [100000]
    if n_defchange > 1:
        orbit_end = np.append(np.delete(orbit_start,0,0), 100000)
    if n_defchange is 0:
        print("oh noes!")
        return(dict)
    for i_change in range(n_defchange):
        if orbit_start[i_change] > orbit or orbit_end[i_change] < orbit:
            continue
        cluspets = cluspets[i_change]
        cluscoad = cluscoad[i_change]
        for i_clus in range(3): #ch8 is last 3 clusters.. always 40 clusters
            clus_start = clusoff1[i_clus]
            clus_end   = clusoff1[i_clus+1]
            pet[clus_start:clus_end] = cluspets[i_clus+37]
            coadd[clus_start:clus_end] = cluscoad[i_clus+37]
    dict['pet'] = pet - petcorr 
    dict['coadd'] = coadd

    fid.close()
    return(dict)

# reads extracted state executions from database (CH8 only)
def read_extracted_states_(orbitrange, state_id, calib_db, in_orbitlist=None, readoutMean=False, readoutNoise=False, orbitList=False, readoutCount=False, errorInTheMean=True):

    if in_orbitlist is None:
        if len(orbitrange) is not 2:
            raise Exception('orbitrange should have 2 elements')

    dict = {}

    dict["orbit_range"] = orbitrange

    #
    # obtain indices to entries of the requested orbit and state
    #

    fid = h5py.File(calib_db, 'r') # generates an exception
    gid = fid["State_"+str('%02d'%state_id)] # generates an exception
    mt_did = gid['metaTable']
    orbitlist = (gid['orbitList'])[:]
    if in_orbitlist is not None:
        #metaindx = orbitlist = in_orbitlist
        metaindx = np.in1d(orbitlist, in_orbitlist)
    else:
        #print('orbitrange=', orbitrange)
        metaindx = (orbitlist >= orbitrange[0]) & (orbitlist <= orbitrange[1])
        #print(metaindx)

    if np.sum(metaindx) == 0:
        raise Exception('orbit range not present in database')

    mtbl = mt_did[metaindx]
    dict['mtbl'] = mtbl

    if readoutMean:
        ds_did      = gid['readoutMean']
        dict['readoutMean'] = ds_did[metaindx,7*1024:]

    if readoutCount:
        ds_did = gid['readoutCount']
        dict['readoutCount'] = ds_did[metaindx,7*1024:]

    # export error-in-the-mean instead of plain stddev.
    if readoutNoise:
        ds_did       = gid['readoutNoise']
        readoutNoise = ds_did[metaindx,7*1024:]
        if errorInTheMean:
            ds_did       = gid['readoutCount']
            readoutNoise /= np.sqrt(ds_did[metaindx,7*1024:])
        readoutNoise *= 1.4826 # sigma = MAD * K 
        dict['readoutNoise'] = readoutNoise

    if orbitList:
        ds_did = gid['orbitList']
        dict['orbitList'] = ds_did[metaindx]

    #
    # find the pet and coadd (TODO: do for multiple orbits)
    #

    orbit = orbitrange[0]
    clusoff1 = [0,10,1014,1024]
    ds_did = gid['clusConf']
    clusconf = ds_did[:]
    orbit_start = clusconf['orbit'][:]
    cluspets = clusconf['pet'][:]
    cluscoad = clusconf['coaddf'][:]
    n_defchange = orbit_start.size 
    pet = np.empty(1024)
    coadd = np.empty(1024)
    #print("CLUSPETS=", cluspets.shape, cluspets)
    if n_defchange is 1:
        orbit_end = [100000]
    if n_defchange > 1:
        orbit_end = np.append(orbit_start[1:], 100000)
    if n_defchange is 0:
        print("oh noes!")
        return(dict)
    for i_change in range(n_defchange):
        if orbit_start[i_change] > orbit or orbit_end[i_change] < orbit:
            continue
        cluspets_ = cluspets[i_change,:]
        cluscoad_ = cluscoad[i_change,:]
        for i_clus in range(3): #ch8 is last 3 clusters.. always 40 clusters
            clus_start = clusoff1[i_clus]
            clus_end   = clusoff1[i_clus+1]
            #print(cluspets_.shape, cluspets_, orbit)
            pet[clus_start:clus_end] = cluspets_[i_clus+37]
            coadd[clus_start:clus_end] = cluscoad_[i_clus+37]
    dict['pet'] = pet - petcorr
    dict['coadd'] = coadd

    fid.close()
    return(dict)

def get_closest_state_exec(orbit, stateid, calib_db, **kwargs):
    fid = h5py.File(calib_db, 'r') # generates an exception
    gid = fid["State_"+str('%02d'%stateid)] # generates an exception
    mt_did = gid['metaTable']
    orbitlist = (gid['orbitList'])[:]
    idx = np.argmin(np.abs(orbitlist[:] - orbit))
    orbit_ = orbitlist[idx]
    print("FOUND:", orbit_)
    fid.close()
    return read_extracted_states_([orbit_,orbit_], stateid, calib_db, **kwargs)

class NoiseModel:
    """
    Class that computed detector pixel noise for specified settings.
    This is a class because various arrays need to be initialized before-hand.
    """

    def __init__(self):
        self.Echarge     = 1.6021892e-19
        self.Boltzman    = 1.38065800E-23
        self.ElBu        = np.array([490.   ,855.   ,854.   ,868.   ,863.   ,529. ,177. ,177. ])               # electrons per BU
        self.Temperature = np.array([200    ,200    ,220    ,220    ,220    ,200  ,150  ,150  ])               # nominal [Kelvin]
        self.Gain        = np.array([40     ,20     ,20     ,20     ,20     ,0    ,0    ,0    ])               # pre amplifier gain
        self.Capacitor   = np.array([4E-12  ,4E-12  ,4E-12  ,4E-12  ,4E-12  ,0    ,0    ,0    ])               # only for ch1-5, [Fahrad]
        self.AdcNoise    = np.zeros(8) + np.sqrt(1./12.)                                                       # bu
        self.AdcThermal  = np.array([0.5    ,0.5    ,0.5    ,0.5    ,0.5    ,0.5  ,0.5  ,0.5  ])               # bu
        self.Opamp       = np.array([500.   ,500.   ,500.    ,500.   ,500.   ,1100., 600. ,600. ]) / self.ElBu # opamp / mux [bu]
        fid = h5py.File("resistance.h5")
        self.resistance = fid["resistance"][:]
        return

    def compute(self, pixel, pet, adc, coaddf=1, give_shotnoise=False):
        """
        Returns noise for specified detector pixel with exposure time pet, and detector filling (excluding analog offset) adc.
        This is computed using on-ground measured resistances.

        Parameters
        ----------

        pixel: int
            pixel number [0..8191]
        pet: float
            pixel exposure time
        adc: float
            detectore filling, excluding analog offset
        coaddf: int, optional
            co-adding factor, default = 1
        give_shotnoise: bool, optional
            if set, also return shotnoise estimate
        """

        channel     = pixel // 1024
        r_junction  = abs(self.resistance[pixel])
        jnoise      = np.sqrt(self.Temperature[channel]*2*self.Boltzman*pet*coaddf/r_junction)/self.Echarge/self.ElBu[channel]           # bu
        ktc         = np.sqrt(self.Boltzman*self.Temperature[channel]*self.Capacitor[channel])/self.Echarge/self.ElBu[channel]                                  # bu
        shotnoise   = np.sqrt(adc/self.ElBu[channel])                                                                                  # bu
        grand_total = np.sqrt(jnoise**2+ktc**2+shotnoise**2+self.AdcNoise[channel]**2+self.AdcThermal[channel]**2+self.Opamp[channel]**2)
        if give_shotnoise:
            return grand_total, shotnoise
        else:
            return grand_total

#-- main -----------------------------------------------------------------------

# test function. not a unit test.. yet
if __name__ == '__main__':
    print("\norbit filter test")
    filt = orbitfilter()
    ra = np.arange(51000)
    mo_mask = filt.get_monthly_orbit_filter(ra)
    qu_mask = filt.get_quality_orbit_filter(ra)
    ra2 = ra[mo_mask]
    ra3 = ra[qu_mask]
    print(ra2.shape, ra3.shape)
    print(filt.get_next_monthly(24000))
    print(filt.get_next_monthly(24044))

    print("\nmemory correction test")
    mc = MemCorrector()
    a = np.arange(8192)*4
    b = mc.correct(a)
    print(b[0:1024])
    print(b[1024:2*1024])
    print(b[2*1024:3*1024])
    print(b[3*1024:4*1024])
    print(b[4*1024:5*1024])

    print("\nnon-linearity correction test")
    nlc = NonlinCorrector()
    b= nlc.correct(a)
    print(b[0:1024])
    print(b[1024:2*1024])
    print(b[2*1024:3*1024])
    print(b[3*1024:4*1024])
    print(b[4*1024:5*1024])
    print(b[5*1024:6*1024])
    print(b[6*1024:7*1024])
    print(b[7*1024:8*1024])

    nm = NoiseModel()
    noise, shot = nm.compute(7*1024+310, 1.0, 5000., give_shotnoise=True)
    pixarr = np.array([7*1024+310,7*1024+311,7*1024+312])
    petarr = np.array([1.0,1.0,1.0])
    adc = np.array([5000.,4500.,4000.])
    noise, shot = nm.compute(pixarr, petarr, adc, give_shotnoise=True)
    print("noise", noise, "shot", shot)
