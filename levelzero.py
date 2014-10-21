from __future__ import print_function, division

import numpy as np
from ctypes import *

#-- globals --------------------------------------------------------------------

ENVI_FILENAME_SIZE = 63
MAX_UTC_STRING = 28

"""
struct mph_envi
{
     char  product[ENVI_FILENAME_SIZE];
     char  proc_stage[2];
     char  ref_doc[24];

     char  acquis[21];
     char  proc_center[7];
     char  proc_time[MAX_UTC_STRING];
     char  soft_version[15];

     char  sensing_start[MAX_UTC_STRING];
     char  sensing_stop[MAX_UTC_STRING];

     char  phase[2];
     short cycle;
     int   rel_orbit;
     int   abs_orbit;
     char  state_vector[MAX_UTC_STRING];
     double delta_ut;
     double x_position;
     double y_position;
     double z_position;
     double x_velocity;
     double y_velocity;
     double z_velocity;
     char  vector_source[3];

     char  utc_sbt_time[MAX_UTC_STRING];
     unsigned int sat_binary_time;
     unsigned int clock_step;

     char  leap_utc[MAX_UTC_STRING];
     short leap_sign;
     char  leap_err[2];

     char  product_err[2];
     unsigned int tot_size;
     unsigned int sph_size;
     unsigned int num_dsd;
     unsigned int dsd_size;
     unsigned int num_data_sets;
};

struct dsd_envi
{
     char name[29];
     char type[2];
     char flname[ENVI_FILENAME_SIZE];
     unsigned int offset;
     unsigned int size;
     unsigned int num_dsr;
     int dsr_size;
};
"""

#-- structures -----------------------------------------------------------------

class dsd_envi(Structure):
	_fields_ = [("name",create_string_buffer(29)),
		("type",create_string_buffer(2)),
		("flname",create_string_buffer(ENVI_FILENAME_SIZE)),
		("offset",c_uint),
		("size",c_uint),
		("num_dsr",c_uint),
		("dsr_size",c_int)]

class mph_envi(Structure):
	_fields_ = [("product",create_string_buffer(ENVI_FILENAME_SIZE)), 
		("proc_stage",create_string_buffer(2)),
		("ref_doc",create_string_buffer(24)),
		("acquis",create_string_buffer(21)),
		("proc_center",create_string_buffer(7)),
		("proc_time",create_string_buffer(MAX_UTC_STRING)),
		("soft_version",create_string_buffer(15)),
		("sensing_start",create_string_buffer(MAX_UTC_STRING)),
		("sensing_stop",create_string_buffer(MAX_UTC_STRING)),
		("phase",create_string_buffer(2)),
		("cycle",c_short),
		("rel_orbit",c_int),
		("abs_orbit",c_int),
		("state_vector",create_string_buffer(MAX_UTC_STRING)),
		("delta_ut",c_double),
		("x_position",c_double),
		("y_position",c_double),
		("z_position",c_double),
		("x_velocity",c_double),
		("y_velocity",c_double),
		("z_velocity",c_double),
	    ("vector_source",create_string_buffer(3)),
	    ("utc_sbt_time",create_string_buffer(MAX_UTC_STRING)),
	    ("sat_binary_time",c_uint),
	    ("clock_step",c_uint),
	    ("leap_utc",create_string_buffer(MAX_UTC_STRING)),
	    ("leap_sign",c_short),
	    ("leap_err",create_string_buffer(2)),
	    ("product_err",create_string_buffer(2)),
	    ("tot_size",c_uint),
	    ("sph_size",c_uint),
	    ("num_dsd",c_uint),
	    ("dsd_size",c_uint),
	    ("num_data_sets",c_uint)]

#-- functions ------------------------------------------------------------------

class LevelZeroFile:

    def __init__(self):
        self.lib = ct.CDLL('levelzerolib.so')

	def read_ch8():

		return 


#-- main -----------------------------------------------------------------------

