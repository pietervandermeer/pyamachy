from __future__ import print_function, division

import numpy as np
import h5py
from ctypes import *
#import numpy.ctypeslib as npct
from envisat import PhaseConverter
from sciamachy_module import NonlinCorrector

#-- globals --------------------------------------------------------------------

ENVI_FILENAME_SIZE = 63
MAX_UTC_STRING = 28
MAX_CLUSTER = 64
NUM_LV0_AUX_BCP = 16
NUM_LV0_AUX_PMTC_FRAME = 5
NUM_LV0_PMD_PACKET = 200

#-- structures -----------------------------------------------------------------

class dsd_envi(Structure):
	_fields_ = [("name",c_char * 29),
		("type",c_char * 2),
		("flname",c_char * ENVI_FILENAME_SIZE),
		("offset",c_uint),
		("size",c_uint),
		("num_dsr",c_uint),
		("dsr_size",c_int)]

class mph_envi(Structure):
	_fields_ = [("product",c_char * ENVI_FILENAME_SIZE),
		("proc_stage",c_char * 2),
		("ref_doc",c_char * 24),
		("acquis",c_char *21),
		("proc_center",c_char * 7),
		("proc_time",c_char * MAX_UTC_STRING),
		("soft_version",c_char * 15),
		("sensing_start",c_char * MAX_UTC_STRING),
		("sensing_stop",c_char * MAX_UTC_STRING),
		("phase",c_char * 2),
		("cycle",c_short),
		("rel_orbit",c_int),
		("abs_orbit",c_int),
		("state_vector",c_char * MAX_UTC_STRING),
		("delta_ut",c_double),
		("x_position",c_double),
		("y_position",c_double),
		("z_position",c_double),
		("x_velocity",c_double),
		("y_velocity",c_double),
		("z_velocity",c_double),
	    ("vector_source",c_char * 3),
	    ("utc_sbt_time",c_char * MAX_UTC_STRING),
	    ("sat_binary_time",c_uint),
	    ("clock_step",c_uint),
	    ("leap_utc",c_char * MAX_UTC_STRING),
	    ("leap_sign",c_short),
	    ("leap_err",c_char * 2),
	    ("product_err",c_char * 2),
	    ("tot_size",c_uint),
	    ("sph_size",c_uint),
	    ("num_dsd",c_uint),
	    ("dsd_size",c_uint),
	    ("num_data_sets",c_uint)]

class sph0_scia(Structure):
	_fields_ = [("descriptor", c_char * 29),
		("start_lat", c_double),
		("start_lon", c_double),
		("stop_lat", c_double),
		("stop_lon", c_double),
		("sat_track", c_double),
		("isp_errors", c_ushort),
		("missing_isps", c_ushort),
		("isp_discard", c_ushort),
		("rs_sign", c_ushort),
		("num_error_isps", c_int),
		("error_isps_thres", c_double),
		("num_miss_isps", c_int),
		("miss_isps_thres", c_double),
		("num_discard_isps", c_int),
		("discard_isps_thres", c_double),
		("num_rs_isps", c_int),
		("rs_thres", c_double),
		("tx_rx_polar", c_char * 6),
		("swath", c_char * 4)]
"""
struct info_clus
{
     unsigned char  chanID;
     unsigned char  clusID;
     unsigned char  coAdding;
     unsigned short start;
     unsigned short length;
};

struct mds0_info
{
     unsigned char   packetID;
     unsigned char   category;
     unsigned char   stateID;
     unsigned char   numClusters;
     unsigned short  length;
     unsigned short  bcps;
     unsigned short  stateIndex;
     unsigned int    offset;
     struct mjd_envi mjd;
     struct info_clus cluster[MAX_CLUSTER];
};

struct mds0_aux
{
     struct mjd_envi   isp;
     struct fep_hdr    fep_hdr;
     struct packet_hdr packet_hdr;
     struct data_hdr   data_hdr;
     struct pmtc_hdr   pmtc_hdr;
     struct aux_src    data_src;
};

struct mds0_det
{
     unsigned short    bcps;
     unsigned short    num_chan;
     int               orbit_vector[8];
     struct mjd_envi   isp;
     struct fep_hdr    fep_hdr;
     struct packet_hdr packet_hdr;
     struct data_hdr   data_hdr;
     struct pmtc_hdr   pmtc_hdr;
     struct det_src    *data_src;
};

struct mjd_envi
{
     int days;
     unsigned int secnd;
     unsigned int musec;
};
"""
class mjd_envi(Structure):
	_fields_ = [("days", c_int),
		("secnd", c_uint),
		("musec", c_uint)]

class info_clus(Structure):
	_fields_ = [("chanID", c_byte),
	    ("clusID", c_byte),
	    ("coAdding", c_byte),
	    ("start", c_ushort),
	    ("length", c_ushort)]

class mds0_info(Structure):
	_fields_ = [("packetID", c_byte),
		("category", c_byte),
		("stateID", c_byte),
		("numClusters", c_byte),
		("length", c_ushort),
		("bcps", c_ushort),
		("stateIndex", c_ushort),
		("offset", c_uint),
		("mjd", mjd_envi),
		("info_clus", info_clus * MAX_CLUSTER)]

"""
struct fep_hdr
{
     struct mjd_envi  gsrt;
     unsigned short   isp_length;
     unsigned short   crc_errs;
     unsigned short   rs_errs;
};
"""

class fep_hdr(Structure):
	_fields_ = [("gsrt", mjd_envi),
		("isp_length", c_ushort),
		("crc_errs", c_ushort),
		("rs_errs", c_ushort)]

"""
struct packet_hdr
{
     union 
     {
	  struct packet_id_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned short op_mode:5;
	       unsigned short vcid:5;
	       unsigned short :6;
#else
	       unsigned short :6;
	       unsigned short vcid:5;
	       unsigned short op_mode:5;
#endif
	  } field;

	  unsigned short two_byte;
     } api;

     unsigned short seq_cntrl;
     unsigned short length;     
};
"""
#puff.. if we don't need these in python.. please don't bother with ctypes unions
#class packet_id_breakout(Union):

class packet_hdr(Structure):
	_fields_ = [("two_byte", c_ushort),
		("seq_cntrl", c_ushort),
		("length", c_ushort)]

"""
struct pmtc_hdr
{
     union 
     {
	  struct pmtc_1_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char ndfm:2;
	       unsigned char :2;
	       unsigned char phase:4;
	       unsigned char sls:2;
	       unsigned char wls:2;
	       unsigned char apsm:2;
	       unsigned char ncwm:2;
#else
	       unsigned char phase:4;
	       unsigned char :2;
	       unsigned char ndfm:2;
	       unsigned char ncwm:2;
	       unsigned char apsm:2;
	       unsigned char wls:2;
	       unsigned char sls:2;
#endif
	  } field;

	  unsigned short two_byte;
     } pmtc_1;
     unsigned short scanner_mode;

     union 
     {
	  struct az_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned short repeat:12;
	       unsigned char  basic:4;
	       unsigned char  h_w:4;
	       unsigned char  rel:4;
	       unsigned char  corr:4;
	       unsigned char  invert:1;
	       unsigned char  filter:1;
	       unsigned char  centre:1;
	       unsigned char  type:1;
#else
	       unsigned char  type:1;
	       unsigned char  centre:1;
	       unsigned char  filter:1;
	       unsigned char  invert:1;
	       unsigned char  corr:4;
	       unsigned char  rel:4;
	       unsigned char  h_w:4;
	       unsigned char  basic:4;
	       unsigned short repeat:12;
#endif
	  } field;

	  unsigned int four_byte;
     } az_param;

     union 
     {
	  struct elv_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned short repeat:12;
	       unsigned char  basic:4;
	       unsigned char  :4;
	       unsigned char  rel:4;
	       unsigned char  corr:4;
	       unsigned char  invert:1;
	       unsigned char  filter:1;
	       unsigned char  centre:1;
	       unsigned char  :1;
#else
	       unsigned char  :1;
	       unsigned char  centre:1;
	       unsigned char  filter:1;
	       unsigned char  invert:1;
	       unsigned char  corr:4;
	       unsigned char  rel:4;
	       unsigned char  :4;
	       unsigned char  basic:4;
	       unsigned short repeat:12;
#endif
	  } field;

	  unsigned int four_byte;
     } elv_param;
     unsigned char factor[6];
};
"""

class pmtc_hdr(Structure):
	_fields_ = [("two_byte1", c_ushort),
		("scanner_mode", c_ushort),
		("four_byte1", c_uint), # 32 bit.. on a 64 bit build? this is bug.. but well.. also note that c_ulong is 64 bit here!
		("four_byte2", c_uint),
		("factor", c_byte * 6)]

"""
struct data_hdr
{
     unsigned char  category;
     unsigned char  state_id;
     unsigned short length;
     union
     {
	  struct vector_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char atc_id:6;
	       unsigned char hsm:2;
	       unsigned char config_id:8;
#else
	       unsigned char hsm:2;
	       unsigned char atc_id:6;
	       unsigned char config_id:8;
#endif
	  } field;

	  unsigned short two_byte;
     } rdv;
     union
     {
	  struct id_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char :4;
	       unsigned char  packet:4;
	       unsigned char  overflow:4;
	       unsigned char :4;
#else
	       unsigned char  packet:4;
	       unsigned char :4;
	       unsigned char :4;
	       unsigned char  overflow:4;
#endif
	  } field;
	  unsigned short two_byte;
     } id;
     unsigned int on_board_time;
};
"""

"""
struct aux_bcp
{
     unsigned short sync;
     unsigned short bcps;
     union
     {
	  struct flag_breakout
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char phase:4;
	       unsigned char m:1;
	       unsigned char d:1;
	       unsigned char eu:1;
	       unsigned char au:1;
	       unsigned char pointing:6;
	       unsigned char :2;
#else
	       unsigned char au:1;
	       unsigned char eu:1;
	       unsigned char d:1;
	       unsigned char m:1;
	       unsigned char phase:4;
	       unsigned char :2;
	       unsigned char pointing:6;
#endif
	  } field;
	  unsigned short two_byte;
     } flags;
     unsigned int azi_encode_cntr;
     unsigned int ele_encode_cntr;
     unsigned short azi_cntr_error;
     unsigned short ele_cntr_error;
     unsigned short azi_scan_error;
     unsigned short ele_scan_error;
};
"""

class aux_bcp(Structure):
	_fields_ = [("sync", c_ushort),
		("bcps", c_ushort),
		("two_byte", c_ushort),
		("azi_encode_cntr", c_uint),
		("ele_encode_cntr", c_uint),
		("azi_cntr_error", c_ushort),
		("ele_cntr_error", c_ushort),
		("azi_scan_error", c_ushort),
		("ele_scan_error", c_ushort)]

class data_hdr(Structure):
	_fields_ = [("category", c_byte),
		("state_id", c_byte),
		("length", c_ushort),
		("two_byte1", c_ushort),
		("two_byte2", c_ushort),
		("on_board_time", c_uint)] # 32 bit, contrary to what you'd expect, most likely a ctypes bug 
#		("on_board_time", c_ulong)]

class pmtc_frame(Structure):
	_fields_ = [("bcp", aux_bcp * NUM_LV0_AUX_BCP),
		("bench_rad", c_ushort), # bench_cntrl
		("bench_elv", c_ushort), # bench_cntrl
		("bench_az", c_ushort)] # bench_cntrl

class mds0_aux(Structure):
	_fields_ = [("isp", mjd_envi),
		("fep_hdr", fep_hdr),
		("packet_hdr", packet_hdr),
		("data_hdr", data_hdr),
		("pmtc_hdr", pmtc_hdr),
		("data_src", pmtc_frame * NUM_LV0_AUX_PMTC_FRAME)] # aux_src structure basically this array of pmtc_frames

"""
struct chan_hdr
{
     unsigned short sync;
     unsigned short bcps;
     unsigned short bias;
     unsigned short temp;
     union
     {
	  struct
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char lu:2;
	       unsigned char is:2;
	       unsigned char id:4;
	       unsigned char clusters:8;
#else
	       unsigned char id:4;
	       unsigned char is:2;
	       unsigned char lu:2;
	       unsigned char clusters:8;
#endif
	  } field;
	  unsigned short two_byte;
     } channel;
     union
     {
	  struct
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char status:3;
	       unsigned char ratio:5;
	       unsigned char frame:8;
#else
	       unsigned char ratio:5;
	       unsigned char status:3;
	       unsigned char frame:8;
#endif
	  } field;
	  unsigned short two_byte;
     } ratio_hdr;
     union
     {
	  struct
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char  cntrl:2;
	       unsigned char  ratio:5;
	       unsigned short sec:9;
	       unsigned char  mode:2;
	       unsigned short etf:14;
#else
	       unsigned short etf:14;
	       unsigned char  mode:2;
	       unsigned short sec:9;
	       unsigned char  ratio:5;
	       unsigned char  cntrl:2;
#endif
	  } field;
	  unsigned int four_byte;
     } command_vis;
     union
     {
	  struct
	  {
#ifdef _SWAP_TO_LITTLE_ENDIAN
	       unsigned char  cntrl:2;
	       unsigned char  pet:4;
	       unsigned char  :2;
	       unsigned char  bias:3;
	       unsigned char  :3;
	       unsigned char  comp:2;
	       unsigned char  mode:2;
	       unsigned short etf:14;
#else
	       unsigned short etf:14;
	       unsigned char  mode:2;
	       unsigned char  comp:2;
	       unsigned char  :3;
	       unsigned char  bias:3;
	       unsigned char  :2;
	       unsigned char  pet:4;
	       unsigned char  cntrl:2;
#endif
	  } field;
	  unsigned int four_byte;
     } command_ir;
};
"""
class chan_hdr(Structure):
	_fields_ = [("sync", c_ushort),
		("bcps", c_ushort),
		("bias", c_ushort),
		("temp", c_ushort),
		("two_byte1", c_ushort),
		("two_byte2", c_ushort),
		("four_byte1", c_uint),
		("four_byte2", c_uint)]
#		("four_byte1", c_ulong),
#		("four_byte2", c_ulong)]

class chan_src(Structure):
	_fields_ = [
		("cluster_id", c_byte),
		("co_adding", c_byte),
		("sync", c_ushort),
		("block_nr", c_ushort),
		("start", c_ushort),
		("length", c_ushort),
		("data", c_char_p)]

class det_src(Structure):
	_fields_ = [("hdr", chan_hdr),
		("pixel", POINTER(chan_src))]

class mds0_det(Structure):
	_fields_ = [("bcps", c_ushort),
		("num_chan", c_ushort),
		("orbit_vector", c_int * 8),
		("isp", mjd_envi),
		("fep_hdr", fep_hdr),
		("packet_hdr", packet_hdr),
		("data_hdr", data_hdr),
		("pmtc_hdr", pmtc_hdr),
		("data_src", POINTER(det_src))]

#-- functions ------------------------------------------------------------------

# convert ctypes mjd_envi to float mjd
def convert_mjd_struct_to_float(struct):
	d = getattr(struct, "days")
	s = getattr(struct, "secnd")
	mu = getattr(struct, "musec")
	return d + (s + mu/1000000.)/(3600*24.)

class LevelZeroFile:
	"""
	class that can read data from scia level 0 files
	"""
	def __init__(self):
		self.lib = CDLL('levelzerolib.so')
		self.mph = mph_envi()
		self.sph = sph0_scia()
		self.pc = PhaseConverter()
		self.nlc = NonlinCorrector()
		self.corrected = False
		self.fname = None 

	def set_file(self, fname):
		self.fname = fname

	# TODO: can better be done with exact char positions, since this is the way l0 filenames are defined
	def get_orbit(self):
		if self.fname is None:
			raise Exception("no filename set!")
		strs = self.fname.split("_")
		if len(strs) < 9:
			raise Exception("wrong amount of underscores in level 0 file!")
		return int(strs[-2])

	def print_mph(self):
		for field_name, field_type in self.mph._fields_:
			print(field_name, getattr(self.mph, field_name))

	def print_sph(self):
		for field_name, field_type in self.sph._fields_:
			print(field_name, getattr(self.sph, field_name))

	def print_dsd(self):
		for i_dsd in range(len(self.dsd)):
			print()
			for field_name, field_type in self.dsd[i_dsd]._fields_:
				print(field_name, getattr(self.dsd[i_dsd], field_name))

	# find the DSD that describes the scia source packets
	def find_scia_source_dsd(self):
		for dsd in self.dsd:
			if getattr(dsd, "name") == "SCIAMACHY_SOURCE_PACKETS":
				return dsd
		return None

	# pre: self.info is present 
	def find_infos(self, packet_id, state_ids):
		infos = []
		for inf in self.info:
			sid = getattr(inf, "stateID")
			if sid in state_ids:
				pid = getattr(inf, "packetID")
				if pid == packet_id:
					infos.append(inf)
		return infos

	def print_info(self, idx):
		for field_name, field_type in self.info[idx]._fields_:
			print(field_name, getattr(self.info[idx], field_name))

	def print_aux(self, idx):
		for field_name, field_type in self.dark_aux[idx]._fields_:
			print(field_name, getattr(self.dark_aux[idx], field_name))

	def print_det(self, idx):
		for field_name, field_type in self.dark_det[idx]._fields_:
			print(field_name, getattr(self.dark_det[idx], field_name))
		orbit_vector = getattr(self.dark_det[idx], "orbit_vector")
		print("orbit_vector:")
		for el in orbit_vector:
			print(el)
		mjd = getattr(self.dark_det[idx], "isp")
		print("isp:")
		for field_name, field_type in mjd._fields_:
			print(field_name, getattr(mjd, field_name))
		fep_hdr = getattr(self.dark_det[idx], "fep_hdr")
		print("fep_hdr:")
		for field_name, field_type in fep_hdr._fields_:
			print(field_name, getattr(fep_hdr, field_name))
		gsrt = getattr(fep_hdr, "gsrt") # fep_hdr.gsrt mjd_envi struct
		for field_name, field_type in gsrt._fields_:
			print(field_name, getattr(gsrt, field_name))
		packet_hdr = getattr(self.dark_det[idx], "packet_hdr")
		print("packet_hdr:")
		for field_name, field_type in packet_hdr._fields_:
			print(field_name, getattr(packet_hdr, field_name))
		data_hdr = getattr(self.dark_det[idx], "data_hdr")
		print("data_hdr:")
		for field_name, field_type in data_hdr._fields_:
			print(field_name, getattr(data_hdr, field_name))
		pmtc_hdr = getattr(self.dark_det[idx], "pmtc_hdr")
		print("pmtc_hdr:")
		for field_name, field_type in pmtc_hdr._fields_:
			print(field_name, getattr(pmtc_hdr, field_name))

	def get_coadd(self, info):
		#info_clus = getattr(info, "info_clus")
		info_clus = info.info_clus
		num_clus = info.numClusters
		print(num_clus)
		return getattr(info_clus[38], "coAdding")

	# output:
	# self.readouts: ch8 detector readouts (1024xN float32)
	# self.dark_det_mjds: MJD2000 (float64) of readouts
	def read_ch8(self, verbose=False):

		#
		# open
		#

		ret = self.lib.OpenFile(self.fname)
		if ret < 0:
 			raise Exception("couldn't open file ["+self.fname+"]. return code = "+str(ret))

		#
		# read main product header
		#

		ret = self.lib._ENVI_RD_MPH(byref(self.mph))
		if ret < 0:
			raise Exception("couldn't read mph")
		if verbose:
			self.print_mph()

		#
		# read special product header
		#

		ret = self.lib._SCIA_LV0_RD_SPH(byref(self.mph), byref(self.sph))
		if ret < 0:
			raise Exception("couldn't read sph")
		if verbose:
			self.print_sph()

		#
		# read dsd's
		#

		n_dsd = self.mph.num_dsd
		if verbose:
			print("allocating "+str(n_dsd)+" dsd records..")
		dsd_arr_type = dsd_envi * n_dsd
		self.dsd = dsd_arr_type()
		if verbose:
			print(self.dsd)
		ret = self.lib._ENVI_RD_DSD(byref(self.mph), byref(self.dsd))
		if ret >= 0:
			n_dsd = ret
		else:
			raise Exception("_ENVI_RD_DSD() failed")
		if verbose:
			self.print_dsd()

		#
		# read all info packets for dark det and source data
		#

		dsd = self.find_scia_source_dsd()
		if dsd is None:
			raise Exception("SCIAMACHY_SOURCE_PACKETS DSD not found!")
		num_dsr = getattr(dsd, "num_dsr")
		info_array_type = mds0_info * num_dsr
		self.info = info_array_type()
		ret = self.lib._SCIA_LV0_RD_MDS_INFO(c_uint(n_dsd), byref(self.dsd), byref(self.info))
		if verbose:
			print("nr info packets =", ret)
		info_dark_det = self.find_infos(1, (8,26,46,63,67)) # 1:DET
		info_dark_aux = self.find_infos(2, (8,26,46,63,67)) # 2:AUX
		n_dark_det = len(info_dark_det)
		n_dark_aux = len(info_dark_aux)
		if verbose:
			print("nr of dark aux =", n_dark_aux)
			print("nr of dark det =", n_dark_det)
		info_dark_det_type = mds0_info * n_dark_det
		info_dark_det_ = info_dark_det_type()
		info_dark_det_[:] = info_dark_det[:]
		info_dark_aux_type = mds0_info * n_dark_aux
		info_dark_aux_ = info_dark_aux_type()
		info_dark_aux_[:] = info_dark_aux[:]

		#
		# read detector data using det info packets
		#

		chan_mask = 1<<(8-1) # channel 8

		det_array_type = mds0_det * n_dark_det
		self.dark_det = det_array_type()

		det_data_type = c_uint * (1024*n_dark_det)
		self.det_data = det_data_type()

		mjds_type = c_double * n_dark_det
		mjds = mjds_type()

		coadds_type = c_ubyte * n_dark_det
		coadds = coadds_type()

		stateids_type = c_ubyte * n_dark_det
		stateids = stateids_type()

		n_readouts = self.lib._SCIA_LV0_RD_DET(byref(info_dark_det_), c_uint(n_dark_det), c_byte(chan_mask), byref(self.dark_det), self.det_data, mjds, coadds, stateids)
		if n_readouts < 0:
			raise Exception("_SCIA_LV0_RD_DET failed! ret="+str(ret))

		#
		# now get rid of useless packets and shrink data buffer
		#

		# turn into 2d numpy array and shrink
		self.readouts = np.frombuffer(self.det_data, dtype=np.uint32).reshape((n_dark_det, 1024))
		#print("TEST", self.readouts[600,0])
		self.readouts = np.array(self.readouts[0:n_readouts, :], dtype=np.float64)
		#print("shape=", self.readouts.shape, "self.readouts=", self.readouts)

		#
		# convert all readout metadata to numpy arrays
		#

		self.mjds = np.array(mjds[0:n_readouts], dtype=np.float64)
		self.stateids = np.array(stateids[0:n_readouts], dtype=np.byte)
		self.coadds = np.array(coadds[0:n_readouts], dtype=np.byte)
		if verbose:
			print("my very own mjds:", self.mjds)

		#
		# convert mjd's to eclipse orbit phase
		#

		phases, orbits = self.pc.get_phase(self.mjds, eclipseMode=True, getOrbits=True)
		self.ephases = orbits + phases
		if verbose:
			print("my very own ephases:", self.ephases)

		#
		# read aux data using aux info packets (really not necessary in this case, but ok.)
		#

		# aux_data_type = mds0_aux * len(info_dark_aux)
		# self.dark_aux = aux_data_type()
		# ret = self.lib._SCIA_LV0_RD_AUX(byref(info_dark_aux_), c_uint(len(info_dark_aux)), self.dark_aux)
		# if ret != len(info_dark_aux):
		# 	raise Exception("_SCIA_LV0_RD_AUX failed!")

		#
		# close
		#

		ret = self.lib.CloseFile()
		if ret < 0:
			raise Exception("couldn't close file. return code = "+str(ret))

	# correct nonlinearity
	def correct_nlin(self):
		self.readouts = self.nlc.correct_ch8(self.readouts)
		self.corrected = True

	# correct for coadding
	def divide_coadd(self):
		n_readouts = self.readouts.shape[0]
#		print(len(self.coadds))
		for i in range(n_readouts):
			self.readouts[i,:] /= self.coadds[i] 
			#print(self.readouts[i,:])

def test_orbit_from_fname(fname):
	l0 = LevelZeroFile()
	test_fname = "/SCIA/LV0_01/O/SCI_NL__0POLRA20040128_161459_000060052023_00427_10000_2049.N1"
	l0.set_file(test_fname)
	print(l0.get_orbit(test_fname))

def saveNoise(orbit, exp_time, noise):
	idx = orbit_dset[:] == orbit
	if idx.sum() > 0:
		# replace
		i = np.nonzero(idx)[0]
		orbit_dset[idx] = orbit
		noise_dset[i,:] = noise
	else:
		# append
		ax1 = orbit_dset.len()
		orbit_dset.resize(ax1+1, axis=0)
		orbit_dset[ax1] = orbit
		ax1 = noise_dset.len()
		noise_dset.resize(ax1+1, axis=0)

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
	import subprocess
	import matplotlib.pyplot as plt
	import warnings
	warnings.simplefilter("error") # warnings to errors
	#np.set_printoptions(threshold=np.nan, precision=4, suppress=True, linewidth=np.nan)

	from vardark_module import load_varkdark_orbit
	from sciamachy_module import get_darkstateid, petcorr, n_chanpix
	import envisat

	n_pix = n_chanpix
	dbname_long = "/SCIA/SDMF31/pieter/vardark_long.h5"

	orbit_range = [10000,10010]
	cmd = "inquire_scia_db.py --orbit="+str(orbit_range[0])+","+str(orbit_range[1])+" --level=0 --proc=O"
#	orbit_range = [44000,44002]
#	cmd = "inquire_scia_db.py --orbit="+str(orbit_range[0])+","+str(orbit_range[1])+" --level=0 --proc=P"

	text_list = subprocess.check_output(cmd, shell=True)
	fnames = text_list.splitlines()

	output_fname = "noise_tng.h5"
	exp_times = [0.125, 0.5, 1.0]
	if h5py.is_hdf5(output_fname):
		fout = h5py.File(output_fname, "r+") # assuming all the groups and datasets are already present.
	else:
		fout = h5py.File(output_fname, "w")
		for exp_time in exp_times:
			print("creating groups..")
			grp = fout.create_group("pet"+str(exp_time))
			orbit_dset = grp.create_dataset("orbit", (0,), maxshape = (envisat.last_orbit,), dtype='i')
			noise_dset = grp.create_dataset("noise", (0, n_pix), maxshape = (envisat.last_orbit, n_pix), dtype='f')

	zero1 = LevelZeroFile()
	# zero.set_file("/SCIA/LV0_01/O/SCI_NL__0POLRA20040128_161459_000060052023_00427_10000_2049.N1")
	for fname in fnames:
		zero2 = LevelZeroFile()
		print(fname)
		zero2.set_file(fname)
		zero2.read_ch8(verbose=False)
		zero2.divide_coadd()
		zero2.correct_nlin() 
		#print(zero2.ephases)
		#print(zero2.dark_det_mjds)
		if zero1.fname is not None:

			#
			# if we have a complete orbit window then get all the readouts that belong to the actual orbit
			#

			orbit = zero1.get_orbit()
			#print("zero1 fname:", zero1.fname)
			#print(zero1.ephases)
			#print(zero1.dark_det_mjds)
			#print(zero1.ephases, zero2.ephases)
			ephases = np.concatenate((zero1.ephases, zero2.ephases))
			readouts = np.concatenate((zero1.readouts, zero2.readouts))
			states = np.concatenate((zero1.stateids, zero2.stateids))
			orbits = np.trunc(ephases) 
			idx = orbits == orbit
			#print(ephases.shape, orbits.shape, idx.shape)
			#print(ephases, orbits, idx)
			ephases = ephases[idx]
			# print(readouts.shape, readouts[600:700,:])
			readouts = readouts[idx,:]
			states = states[idx]
			wave_phases, wave_orbit, wave, wave_ao = load_varkdark_orbit(orbit, False)
			wave_phases = np.mod(wave_phases, 1)

			#
			# compute and store noise for each exposure time (TODO: separate function?)
			#

			for exp_time in exp_times:

				print("exp_time=", exp_time)

				#
				# open output group and its datasets
				#

				grp = fout["pet"+str(exp_time)]
				orbit_dset = grp["orbit"]
				noise_dset = grp["noise"]

				#
				# get all readouts for the state id that matches the exposure time
				#

				stateid = get_darkstateid(exp_time, orbit)
				idx = states == stateid
				ephases_ = np.mod(ephases[idx],1)
				readouts_ = np.array(readouts[idx,:], dtype=np.float64)

				#
				# subtract dark signal to remove any form of trend
				#

				i_phase = 0
				for phase in ephases_:  # TODO phase interpolation
					ephase_idx = np.argmin(np.abs(phase - wave_phases))
					diffs = (phase - wave_phases) * 50
					lephase_idx = np.where((diffs >= 0.) & (diffs <= 1.))[0]
					if lephase_idx.size == 0:
						lephase_idx = ephase_idx
						frac = 0
					elif lephase_idx.size > 1:
						lephase_idx = lephase_idx[0]
						frac = diffs[lephase_idx]
					else:
						frac = diffs[lephase_idx]
					interp_wave = frac*wave[0,lephase_idx+1,:] + (1-frac)*wave[0,lephase_idx,:]
					interp_wave = interp_wave.flatten()
					# 0 order
#					readouts_[i_phase,:] -= wave[ephase_idx,:]*(exp_time-petcorr) + wave_ao.reshape((n_pix))
					# linear interpolation
					readouts_[i_phase,:] -= interp_wave*(exp_time-petcorr) + wave_ao.reshape((n_pix))
					i_phase += 1

				if (exp_time == 1.0):
					pixnr=115
					plt.cla()
					plt.plot(ephases_, readouts_[:,pixnr], 'bo', label='corr readouts')
					plt.plot([0,.25], [0,0], label='0')
					ranges = [[0,.10],[.11,.17],[.18,.24]]
					stds = np.empty((3))
					cents = np.empty((3))
					mds = np.empty((3))
					i = 0
					for rng in ranges:
						idx = (ephases_ > rng[0]) & (ephases_ < rng[1])
						mn = np.mean(readouts_[idx,pixnr])
						md = np.median(readouts_[idx,pixnr])
						mds[i] = md
						stds[i] = np.std(readouts_[idx,pixnr]) / np.sqrt(np.sum(idx))
						cents[i] = np.mean(ephases_[idx])
						#plt.plot([rng[0],rng[1]], [mn,mn], label='mean '+str(rng[0]))
						#plt.plot([rng[0],rng[1]], [md,md], label='median '+str(rng[0]))
						i+=1
					plt.errorbar(cents, mds, stds, fmt='ro', ls='none', label='median,sig-mu')
					plt.legend(loc='best')
					plt.show()
				#print()

				#
				# compute noise figures
				#

				noise = readouts_.std(axis=0)
				#print("NOISE", noise)
				# TODO: median abs dev?
				# TODO: also compute mean to verify there is no bias? 

				#
				# store noise figures
				#

				#print("noise=", noise)
				#print("readouts=", readouts)
				saveNoise(orbit, exp_time, noise)

		# prepare for next orbit
		zero1 = zero2

	fout.close()
