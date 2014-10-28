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
FILE *fd_nadc;
bool Use_Extern_Alloc = true;
bool File_Is_Open = false;

/*+++++ Static Variables +++++*/
static const char err_msg[] = "invalid number of function arguments";

//#define DEBUG

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
       fprintf(stderr, "\n" );
       goto done;
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
	  fprintf(stderr, "No open stream\n" );
      goto done;
     }
     num_info = SCIA_LV0_RD_MDS_INFO(fd_nadc, num_dsd, dsd, &info);
     if ( IS_ERR_STAT_FATAL ) {
	  fprintf(stderr, "GET_SCIA_LV0_MDS_INFO\n" );
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
	  fprintf(stderr, "No open stream\n" );
      goto done;
     }

     nr_aux = (int) SCIA_LV0_RD_AUX( fd_nadc, info, num_info, &aux );
     if ( IS_ERR_STAT_FATAL ) return -1;

     return nr_aux;
 done:
     return -1;
}

int _SCIA_LV0_RD_DET (struct mds0_info *info, unsigned int num_det, unsigned char chan_mask, struct mds0_det *C_det, unsigned int *data, double *mjds, uint8_t *coadds, uint8_t *stateids)
{
    const char prognm[] = "_SCIA_LV0_RD_DET";
    unsigned int n_ch, n_cl;
    unsigned int  num_clus;
    unsigned int nr = 0;
    unsigned int offs = 0;
    unsigned int last_offs = 0;
    unsigned int sz_data = 0;
    unsigned int nrec_out = 0;

    if ( fileno( fd_nadc ) == -1 ) 
    {
        NADC_ERROR( prognm, NADC_ERR_FILE, "No open stream" );
        return -1;
    }
    struct det_src *tmp_buf = malloc(10*sizeof(struct det_src));

    nadc_stat = NADC_STAT_SUCCESS;
    for ( nr = 0; nr < num_det; nr++, C_det++ ) 
    {
#ifdef DEBUG
        printf("nr=%d\n", nr);
#endif

        for ( n_cl = 0; n_cl < (unsigned int) info[nr].numClusters; n_cl++ )
            sz_data += info[nr].cluster[n_cl].length;

        C_det->data_src = tmp_buf;
        SCIA_LV0_RD_DET( fd_nadc, info+nr, 1, chan_mask, &C_det );

        if ( IS_ERR_STAT_FATAL ) 
        {
            fprintf(stderr, "fatal!\n");
            return -1;
        }

#ifdef DEBUG
        printf("nr=%d, nchan=%d\n", nr, C_det->num_chan); 
#endif
        for ( n_ch = 0; n_ch < C_det->num_chan; n_ch++ ) 
        {
            num_clus = C_det->data_src[n_ch].hdr.channel.field.clusters ;
#ifdef DEBUG
            printf("n_ch=%d, num_clus=%d\n", n_ch, num_clus);
#endif
            for ( n_cl = 0; n_cl < num_clus; n_cl++ ) 
            {
#ifdef DEBUG
                printf("offs=%d\n", offs);
#endif
                if ( (offs+C_det->data_src[n_ch].pixel[n_cl].length) > sz_data ) 
                {
                    fprintf(stderr, "size?!\n");
                    goto done;
                }
                UNPACK_LV0_PIXEL_VAL( &C_det->data_src[n_ch].pixel[n_cl], data+offs );
                offs += C_det->data_src[n_ch].pixel[n_cl].length;
            } // cluster loop

        } // channel loop

        if (offs - last_offs == 1024)
        {
            // convert mjd struct to floating point, also add the bcps to get the real time of the readout
            struct mjd_envi isp = C_det->isp;
            uint16_t bcps = C_det->bcps;
            mjds[nrec_out] = isp.days + (isp.secnd + (isp.musec/1000000.) + (bcps/16.))/(3600*24.);

            // get the coadding factor
            int nclusters = info[nr].numClusters;
            coadds[nrec_out] = info[nr].cluster[nclusters-2].coAdding;
            stateids[nrec_out] = info[nr].stateID;

#if DEBUG
            int nn;
            for (nn=0; nn<nclusters; nn++)
                printf("%d ", info[nr].cluster[nn].coAdding);
            printf("\n");
            printf("nrec_out=%d, mjd=%f, coadd=%d, stateid=%d\n", nrec_out, mjds[nrec_out], coadds[nrec_out], stateids[nrec_out]);
#endif

            nrec_out++;
        }

        last_offs = offs;
    } // det packet loop
    free(tmp_buf);
    return nrec_out;
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

// returns all ch8 det readouts and mjd's
int _GET_ALL_CH8_DET(char *fname, float *readouts, int *n_readouts)
{
    unsigned int indx_dsd;
    unsigned int num_info;
    unsigned int num_dsd;
    struct mph_envi  mph;
    struct sph0_scia sph;
    struct mds0_info *info;
    struct dsd_envi *dsd;
    int ret;
    struct mds0_aux *aux_out;
    const char dsd_name[] = "SCIAMACHY_SOURCE_PACKETS";
    FILE *fd;

    if ((fd = fopen( fname, "r" )) == NULL ) {
        fprintf(stderr, strerror( errno ) );
        fprintf(stderr, "\n" );
        goto done;
    }

    // read Main Product Header
    ENVI_RD_MPH( fd, &mph );
    //num_dsd = mph.num_dsd;

    // read special product header
    SCIA_LV0_RD_SPH(fd, mph, &sph);

    num_dsd = (int) ENVI_RD_DSD( fd_nadc, mph, dsd );
    if ( IS_ERR_STAT_FATAL ) {
        fprintf(stderr, "ENVI_RD_DSD() failed!\n");
    }

    // get index to data set descriptor
    indx_dsd = GET_ENVI_DSD_INDX( num_dsd, dsd, dsd_name );
    if ( IS_ERR_STAT_ABSENT ) {
        fprintf(stderr, dsd_name );
        goto done;
    }

    // allocate memory to store info-records
    num_info = dsd[indx_dsd].num_dsr;
    info = (struct mds0_info *) malloc( (size_t) num_info * sizeof(struct mds0_info) );
    if ( info == NULL ) {
        fprintf(stderr, "error mds0_info\n");
        goto done;
    }

    // extract info-records from input file
    num_info = GET_SCIA_LV0_MDS_INFO( fd, mph, &dsd[indx_dsd], info );
    printf("num_info = %d\n", num_info);

    // TODO: filter on state id  (8,26,46,63,67) and packet_id (AUX)

    int num_aux, i_info;
    for (i_info = 0, num_aux = 0; i_info < num_info; i_info++)
    {
        struct mds0_info *inf;
        inf = info+i_info;
        if ((inf->stateID == 8) || (inf->stateID == 26) || (inf->stateID == 46) || (inf->stateID == 63) || (inf->stateID == 63))
        {
            if (inf->packetID == 2) // AUX
            {
                memmove(info+num_aux, info+i_info, sizeof(struct mds0_info));
                num_aux++;
            }
        }
    }

    // read Auxiliary source packets
    SCIA_LV0_RD_AUX(fd, info, num_info, &aux_out);
    //IF status NE 0 THEN NADC_ERR_TRACE

    // read Detector source packets
    //PRINT, 'SCIA_LV0_RD_MDS_DET'
    //SCIA_LV0_RD_DET, info, mds_det, status=status, state_id=6
    //IF status NE 0 THEN NADC_ERR_TRACE

    ret = fclose( fd );
    if (ret < 0) {
        fprintf(stderr, "failed to close file\n" );
        goto done;
    }

    return 0;
done:
    return -1;        
}
