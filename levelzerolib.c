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
bool File_Is_Open = false;

/*+++++ Static Variables +++++*/
static const char err_msg[] = "invalid number of function arguments";

/*+++++++++++++++++++++++++ Static Functions +++++++++++++++++++++++*/
static inline
void UNPACK_LV0_PIXEL_VAL( const struct chan_src *pixel,
                           /*@out@*/ unsigned int *data )
{
     register unsigned short np = 0;

     register unsigned char *cpntr = pixel->data;

     if ( pixel->co_adding == (unsigned char) 1 ) {
       do {
            *data++ = (unsigned int) cpntr[1]
              + ((unsigned int) cpntr[0] << 8);
            cpntr += 2;
       } while ( ++np < pixel->length );
     } else {
       do {
            *data++ = (unsigned int) cpntr[2]
              + ((unsigned int) cpntr[1] << 8)
              + ((unsigned int) cpntr[0] << 16);
            cpntr += 3;
       } while ( ++np < pixel->length );
     }
}

/*+++++++++++++++++++++++++ Main Program or Functions +++++++++++++++*/

int OpenFile( char* fname )
{
     if ( File_Is_Open ) (void) fclose( fd_nadc );

     if ((fd_nadc = fopen( fname, "r" )) == NULL ) {
       fprintf(stderr, strerror( errno ) );
     } else
       File_Is_Open = true;

     return 0;
 done:
     return -1;
}

int CloseFile( void )
{
     int stat = 0;

     if ( File_Is_Open ) {
       stat = fclose( fd_nadc );
       File_Is_Open = FALSE;
     }
     return stat;
}

int _ENVI_RD_MPH( struct mph_envi *mph )
{
     const char prognm[] = "_ENVI_RD_MPH";

     if ( ! File_Is_Open ) {
       fprintf(stderr, "No open stream" );
       goto done;
     }

     ENVI_RD_MPH( fd_nadc, mph );
     if ( IS_ERR_STAT_FATAL ) return -1;

     return 1;
 done:
     return -1;
}

int _ENVI_RD_DSD( struct mph_envi *pmph, struct dsd_envi *dsd )
{
     const char prognm[] = "_ENVI_RD_DSD";

     int nr_dsd = 0;

     struct mph_envi  mph;

     if ( ! File_Is_Open ) {
       fprintf(stderr, "No open stream" );
       goto done;
     }

     mph = *pmph;

     nr_dsd = (int) ENVI_RD_DSD( fd_nadc, mph, dsd );
     if ( IS_ERR_STAT_FATAL ) return -1;

     return nr_dsd;
 done:
     return -1;
}

int _SCIA_LV0_RD_SPH (struct mph_envi *pmph, struct sph0_scia *sph)
{
     const char prognm[] = "_SCIA_LV0_RD_SPH";

     struct mph_envi  mph;

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
	  fprintf(stderr, "No open stream" );
      goto done;
     }

     nr_aux = (int) SCIA_LV0_RD_AUX( fd_nadc, info, num_info, &aux );
     if ( IS_ERR_STAT_FATAL ) return -1;

     return nr_aux;
 done:
     return -1;
}

int _SCIA_LV0_RD_DET (struct mds0_info *info, unsigned int num_det, unsigned char chan_mask, struct mds0_det *C_det, unsigned int *data)
{
     const char prognm[] = "_SCIA_LV0_RD_DET";

     register unsigned int n_ch, n_cl;

     unsigned int  num_clus;

     register unsigned int nr = 0;
     register unsigned int offs = 0;

     unsigned int sz_data = 0;

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

	  for ( n_ch = 0; n_ch < C_det->num_chan; n_ch++ ) {
	       num_clus = C_det->data_src[n_ch].hdr.channel.field.clusters ;

	       for ( n_cl = 0; n_cl < num_clus; n_cl++ ) {

		    if ( (offs+C_det->data_src[n_ch].pixel[n_cl].length)
			 > sz_data ) goto done;
		    UNPACK_LV0_PIXEL_VAL( &C_det->data_src[n_ch].pixel[n_cl],
					  data+offs );
		    offs += C_det->data_src[n_ch].pixel[n_cl].length;
	       }
	  }
     }
     return num_det;
 done:
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
