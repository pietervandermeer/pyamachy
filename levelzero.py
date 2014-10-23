from __future__ import print_function, division

import numpy as np
from ctypes import *
import numpy.ctypeslib as npct

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

class LevelZeroFile:
	"""
	class that can read data from scia level 0 files
	"""
	def __init__(self):
		self.lib = CDLL('levelzerolib.so')
		self.mph = mph_envi()
		self.sph = sph0_scia()

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

	def read_ch8(self):
		#
		# open
		#

		ret = self.lib.OpenFile("/SCIA/LV0_01/O/SCI_NL__0POLRA20040128_161459_000060052023_00427_10000_2049.N1")
		if ret < 0:
			raise Exception("couldn't open file. return code = "+str(ret))

		#
		# read main product header
		#

		ret = self.lib._ENVI_RD_MPH(byref(self.mph))
		if ret < 0:
			raise Exception("couldn't read mph")
		self.print_mph()

		#
		# read special product header
		#

		ret = self.lib._SCIA_LV0_RD_SPH(byref(self.mph), byref(self.sph))
		if ret < 0:
			raise Exception("couldn't read sph")
		self.print_sph()

		#
		# read dsd's
		#

		n_dsd = self.mph.num_dsd
		print("allocating "+str(n_dsd)+" dsd records..")
		dsd_arr_type = dsd_envi * n_dsd
		self.dsd = dsd_arr_type()
		print(self.dsd)
		ret = self.lib._ENVI_RD_DSD(byref(self.mph), byref(self.dsd))
		print(ret)
		if ret >= 0:
			n_dsd = ret
		else:
			raise Exception("_ENVI_RD_DSD() failed")
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
		print("nr info packets =", ret)
		#self.print_info(5102)
		info_dark_det = self.find_infos(1, (8,26,46,63,67)) # 1:DET
		info_dark_aux = self.find_infos(2, (8,26,46,63,67)) # 2:AUX
		print("nr of dark aux =", len(info_dark_aux))
		print("nr of dark det =", len(info_dark_det))
		info_dark_det_type = mds0_info * len(info_dark_det)
		info_dark_det_ = info_dark_det_type()
		info_dark_det_[:] = info_dark_det[:]
		info_dark_aux_type = mds0_info * len(info_dark_aux)
		info_dark_aux_ = info_dark_aux_type()
		info_dark_aux_[:] = info_dark_aux[:]

		#
		# read detector data using det info packets
		#

		chan_mask = 1<<(8-1) # channel 8
		det_array_type = mds0_det * len(info_dark_det)
		self.dark_det = det_array_type()
		det_data_type = c_uint * (1024*len(info_dark_det))
		self.mountain = det_data_type()

		# allocate the det_src structs that are pointed to by mds0_det
		# although.. forget it. too much hassle. and it seems to be temporary in nature as well! 
# 		det_src_type = det_src * len(info_dark_det)
# 		data_src = det_src_type()
# 		i = 0
# 		for d in self.dark_det:
# 			#print(data_src[i])
# #			setattr(d, "data_src", POINTER(det_src)())
# 			data_src = (det_src * 8) ()
# 			print(data_src)
# 			setattr(d, "data_src", POINTER(data_src[i]))
# 			i += 1

		ret = self.lib._SCIA_LV0_RD_DET(byref(info_dark_det_), c_uint(len(info_dark_det)), c_byte(chan_mask), byref(self.dark_det), self.mountain)
		if ret != len(info_dark_det):
			raise Exception("_SCIA_LV0_RD_DET failed! ret="+str(ret))
		print(self.mountain[1024:2048])
		print(len(info_dark_det))
		self.print_det(1000)
		print()

		#
		# read aux data using aux info packets
		#

		aux_data_type = mds0_aux * len(info_dark_aux)
		self.dark_aux = aux_data_type()
		ret = self.lib._SCIA_LV0_RD_AUX(byref(info_dark_aux_), c_uint(len(info_dark_aux)), self.dark_aux)
		if ret != len(info_dark_aux):
			raise Exception("_SCIA_LV0_RD_AUX failed!")
		print()
		self.print_aux(0)
		print()

		#
		# close
		#

		ret = self.lib.CloseFile()
		if ret < 0:
			raise Exception("couldn't close file. return code = "+str(ret))

#-- main -----------------------------------------------------------------------

if __name__ == "__main__":
	zero = LevelZeroFile()
	zero.read_ch8()
