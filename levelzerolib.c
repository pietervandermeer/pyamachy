/* based on code by richard van hees for nadc tools (idl interface).. recoded to easily link to python code */
/*
 * Define _POSIX_SOURCE to indicate
 * that this is a POSIX program
 */
#define  _POSIX_SOURCE 2

/*+++++ System headers +++++*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>

/*+++++ Local Headers +++++*/
#define _SCIA_LEVEL_0
#include <defs_nadc.h>

/*+++++ Global Variables +++++*/
extern FILE *fd_nadc;

/*+++++ Static Variables +++++*/
static const char err_msg[] = "invalid number of function arguments";

/*+++++++++++++++++++++++++ Static Functions +++++++++++++++++++++++*/
#include <../libSCIA/Inline/unpack_lv0_data.inc>

/*+++++++++++++++++++++++++ Main Program or Functions +++++++++++++++*/
int _SCIA_LV0_RD_SPH (struct mph_envi *pmph, struct sph0_scia *sph)
{
     const char prognm[] = "_SCIA_LV0_RD_SPH";

     struct mph_envi  mph;
     struct sph0_scia *sph;

     if ( fileno( fd_nadc ) == -1 ) {
	  fprintf(stderr, "No open stream");
      goto done;
     }
     mph = *pmph;

     SCIA_LV0_RD_SPH( fd_nadc, mph, sph );
     if ( IS_ERR_STAT_FATAL )
        return -1;

     return 1;
 done:
     return -1;
}

int _GET_LV0_MDS_INFO(struct mph_envi *pmph, struct dsd_envi *pdsd, struct mds0_info *info)
{
     const char prognm[] = "_GET_LV0_MDS_INFO";

     unsigned int num_info;

     struct mph_envi  mph;
     struct dsd_envi  dsd;
     struct mds0_info *info;

     if ( fileno( fd_nadc ) == -1 ) {
	  fprintf(stderr, "No open stream" );
      goto done;
     }
     mph = *pmph;
     dsd = *pdsd;

     num_info = GET_SCIA_LV0_MDS_INFO( fd_nadc, mph, &dsd, info );
     if ( IS_ERR_STAT_FATAL ) {
	  fprintf(stderr, "GET_SCIA_LV0_MDS_INFO" );
      goto done;
     }
     return (int) num_info;
 done:
     return -1;
}


int _SCIA_LV0_RD_MDS_INFO (unsigned int num_dsd, struct dsd_envi *dsd, struct mds0_info *info)
{
     const char prognm[] = "_SCIA_LV0_RD_MDS_INFO";

     unsigned int num_info;
     
     struct dsd_envi  *dsd;
     struct mds0_info *info;

     if ( fileno( fd_nadc ) == -1 ) {
	  fprintf(stderr, "No open stream" );
      goto done;
     }
     num_info = SCIA_LV0_RD_MDS_INFO( fd_nadc, num_dsd, dsd, &info );
     if ( IS_ERR_STAT_FATAL ) {
	  fprintf(stderr, "GET_SCIA_LV0_MDS_INFO" );
      goto done;
     }

     return (int) num_info;
 done:
     return -1;
}

int _SCIA_LV0_RD_AUX ( struct mds0_info *info, unsigned int num_info, struct mds0_aux *aux)
{
     const char prognm[] = "_SCIA_LV0_RD_AUX";

     int nr_aux;

     if ( fileno( fd_nadc ) == -1 ) {
	  fprintf(stderr, NADC_ERR_FILE, "No open stream" );
      goto done;
     }

     nr_aux = (int) SCIA_LV0_RD_AUX( fd_nadc, info, num_info, &aux );
     if ( IS_ERR_STAT_FATAL ) return -1;

     return nr_aux;
 done:
     return -1;
}

int _SCIA_LV0_RD_DET (struct mds0_info *info, unsigned int num_det, signed char cbuff, struct IDL_mds0_det *det, unsigned int *data)
{
     const char prognm[] = "_SCIA_LV0_RD_DET";

     register unsigned int n_ch, n_cl;

     unsigned char chan_mask;
     unsigned int  num_clus;

     register unsigned int nr = 0;
     register unsigned int offs = 0;

     unsigned int sz_data = 0;
     unsigned int *data;

     struct mds0_det  *C_det;
     struct mds0_info *info;

     struct IDL_chan_src
     {
	  unsigned char  cluster_id;
	  unsigned char  co_adding;
	  unsigned short sync;
	  unsigned short block_nr;
	  unsigned short start;
	  unsigned short length;
	  IDL_ULONG      pntr_data;            /* IDL uses 32-bit addresses */
     };

     struct IDL_det_src
     {
	  struct chan_hdr hdr;
	  struct IDL_chan_src pixel[LV0_MX_CLUSTERS];

     };

     struct IDL_mds0_det 
     { 
	  unsigned short      bcps;
	  unsigned short      num_chan;
	  int                 orbit_vector[8];
	  struct mjd_envi     isp;
	  struct fep_hdr      fep_hdr;
	  struct packet_hdr   packet_hdr;
	  struct data_hdr     data_hdr;
	  struct pmtc_hdr     pmtc_hdr;
	  struct IDL_det_src  data_src[SCIENCE_CHANNELS];
     } *det;

     if ( fileno( fd_nadc ) == -1 ) {
	  NADC_ERROR( prognm, NADC_ERR_FILE, "No open stream" );
	  return -1;
     }
#if 0
     info  = (struct mds0_info *) argv[0];
     num_det   = *(unsigned int *) argv[1];
     cbuff = *(char *) argv[2];
     chan_mask = (unsigned char) cbuff;
     det   = (struct IDL_mds0_det *) argv[3];
     data  = (unsigned int *) argv[4];
#endif
/*
 * we need SCIA_LV0_RD_DET to do dynamic memory allocation
 */
     nadc_stat = NADC_STAT_SUCCESS;
     for ( nr = 0; nr < num_det; nr++ ) {
	  for ( n_cl = 0; n_cl < (unsigned int) info[nr].numClusters; n_cl++ )
	       sz_data += info[nr].cluster[n_cl].length;
/*
 * read one detector MDS into memory
 */
	  Use_Extern_Alloc = FALSE;
	  SCIA_LV0_RD_DET( fd_nadc, info+nr, 1, chan_mask, &C_det );
	  Use_Extern_Alloc = TRUE;
	  if ( IS_ERR_STAT_FATAL ) return -1;
/*
 * copy C-struct to IDL-struct
 */
	  det[nr].bcps = C_det->bcps;
	  det[nr].num_chan = C_det->num_chan; 
	  (void) memcpy( det[nr].orbit_vector, C_det->orbit_vector, 
			 8 * sizeof( int  ) );
	  (void) memcpy( &det[nr].isp, &C_det->isp, 
			 sizeof( struct mjd_envi ) );
	  (void) memcpy( &det[nr].fep_hdr, &C_det->fep_hdr, 
			 sizeof( struct fep_hdr ) );
	  (void) memcpy( &det[nr].packet_hdr, &C_det->packet_hdr, 
			 sizeof( struct packet_hdr ) );
	  (void) memcpy( &det[nr].data_hdr, &C_det->data_hdr, 
			 sizeof( struct data_hdr ) );
	  (void) memcpy( &det[nr].pmtc_hdr, &C_det->pmtc_hdr, 
			 sizeof( struct pmtc_hdr ) );

	  for ( n_ch = 0; n_ch < C_det->num_chan; n_ch++ ) {
	       num_clus = C_det->data_src[n_ch].hdr.channel.field.clusters ;

	       (void) memcpy( &det[nr].data_src[n_ch].hdr,
                              &C_det->data_src[n_ch].hdr,
			      sizeof( struct chan_hdr ) );

	       for ( n_cl = 0; n_cl < num_clus; n_cl++ ) {
		    det[nr].data_src[n_ch].pixel[n_cl].cluster_id =
			 C_det->data_src[n_ch].pixel[n_cl].cluster_id;
		    det[nr].data_src[n_ch].pixel[n_cl].co_adding =
			 C_det->data_src[n_ch].pixel[n_cl].co_adding;
		    det[nr].data_src[n_ch].pixel[n_cl].sync =
			 C_det->data_src[n_ch].pixel[n_cl].sync;
		    det[nr].data_src[n_ch].pixel[n_cl].block_nr =
			 C_det->data_src[n_ch].pixel[n_cl].block_nr;
		    det[nr].data_src[n_ch].pixel[n_cl].start =
			 C_det->data_src[n_ch].pixel[n_cl].start;
		    det[nr].data_src[n_ch].pixel[n_cl].length =
			 C_det->data_src[n_ch].pixel[n_cl].length;

		    if ( (offs+C_det->data_src[n_ch].pixel[n_cl].length)
			 > sz_data ) goto done;
		    UNPACK_LV0_PIXEL_VAL( &C_det->data_src[n_ch].pixel[n_cl],
					  data+offs );
		    offs += C_det->data_src[n_ch].pixel[n_cl].length;
	       }
	  }
/*
 * release memory
 */
	  SCIA_LV0_FREE_MDS_DET( 1, C_det );
     }
     return num_det;
 done:
     SCIA_LV0_FREE_MDS_DET( 1, C_det );
     return -1;
}

#if 0
int IDL_STDCALL _SCIA_LV0_RD_PMD ( int argc, void *argv[] )
{
     const char prognm[] = "_SCIA_LV0_RD_PMD";

     int num_info, nr_pmd;

     struct mds0_info *info;
     struct mds0_pmd *pmd;

     if ( argc != 3 ) NADC_GOTO_ERROR( prognm, NADC_ERR_PARAM, err_msg );
     if ( fileno( fd_nadc ) == -1 ) 
	  NADC_GOTO_ERROR( prognm, NADC_ERR_FILE, "No open stream" );

     info  = (struct mds0_info *) argv[0];
     num_info  = *(int *) argv[1];
     pmd   = (struct mds0_pmd *) argv[2];

     nr_pmd = (int) SCIA_LV0_RD_PMD( fd_nadc, info, num_info, &pmd );
     if ( IS_ERR_STAT_FATAL ) return -1;

     return nr_pmd;
 done:
     return -1;
}
#endif
