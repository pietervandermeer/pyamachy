#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <defs_nadc.h>
#include <defs_scia.h>
#include <proto_scia.h>

bool Use_Extern_Alloc = true;

void myprint(void);

void print_doubles(double *ds, long n)
{
    long i;
    for (i=0; i<n; i++)
    {
        printf("%f\n", ds[i]);
    }
}

void print_double(double d)
{
    printf("%f\n", d);
}

void myprint()
{
    printf("hello world\n");
}

struct roe_rec {
     double         julianDay;
     unsigned short orbit;
     unsigned short relOrbit;
     unsigned char  phase;
     unsigned char  cycle;
     unsigned short repeat;
     unsigned char  saaDay;
     unsigned char  saaEclipse;
     double         eclipseExit;
     double         eclipseEntry;
     double         period;
     double         anxLongitude;
     char           UTC_anx[28];
     char           UTC_flt[28];
     char           mlst[8];
};

/*+++++++++++++++++++++++++
NOTE: modified version of the original nadc_tools function!! orbitPhase is not made absolute!
.IDENTifer   GET_SCIA_ROE_INFO
.PURPOSE     get orbit phase and SAA flag from ROE records
.INPUT/OUTPUT
  call as   GET_SCIA_ROE_INFO( eclipseMode, jday, 
                               &absOrbit, &saaFlag, &orbitPhase );
     input:
             bool   eclipseMode :  TRUE  - orbit phase used for SDMF (v2.4)
                                   FALSE - orbit phase as used by ESA
             double jday        :  Julian day (# days since 2000-01-01)
    output:
             int   *absOrbit    :  absolute orbit number
             float *saaFlag     :  in-precise SAA flag
             float *orbitPhase  :  orbit phase [0,1]

.RETURNS     nothing
             error status passed by global variable ``nadc_stat''
.COMMENTS    none
-------------------------*/
void GET_SCIA_ROE_INFO_MOD( bool eclipseMode, const double jday, 
                        int *absOrbit, bool *saaFlag, float *orbitPhase )
{
     const double secPerDay = 60. * 60 * 24;

     static struct roe_rec roe;

     double tdiff, phase;
/*
 * initialise return values
 */
     *absOrbit   = -1;
     *saaFlag    = FALSE;
     *orbitPhase = -1.f;
/*
 * get absOrbit and ROE entry for given julianDay
 */
     if ( _GET_ROE_ENTRY( jday, &roe ) < 0 ) return;

     /* set absolute orbit number */
     *absOrbit = roe.orbit;

     /* calculate orbit phase */
     tdiff = (jday - roe.julianDay) * secPerDay - roe.eclipseEntry;
     phase = fmod( tdiff, roe.period ) / roe.period;
     if ( ! eclipseMode ) {
          phase += 0.5 *
               ((roe.eclipseEntry - roe.eclipseExit) / roe.period - 0.5);
     }
     *orbitPhase = phase; // ((phase >= 0) ? (float) phase : (float) (phase + 1));

     /* set SAA flag */
     tdiff = (jday - roe.julianDay) * secPerDay;
     if ( tdiff > roe.eclipseExit && tdiff < roe.eclipseEntry )
          *saaFlag = (roe.saaDay == UCHAR_ZERO);
     else
          *saaFlag = (roe.saaEclipse == UCHAR_ZERO);
}

void _GET_SCIA_ROE_ORBITPHASE(bool eclipseMode, bool saaFlag, int32_t num, int32_t orbit, float *orbitPhase, double *julianDay)
{
    int32_t nm;
    for (nm = 0; nm < num; nm++)
    {
        GET_SCIA_ROE_INFO_MOD(eclipseMode, julianDay[nm],
                              &orbit, &saaFlag, orbitPhase+nm);
    }
}

void _GET_SCIA_ROE_ORBITPHASE_ORBIT(bool eclipseMode, bool saaFlag, int32_t num, int32_t *orbit, float *orbitPhase, double *julianDay)
{
    int32_t nm;
    for (nm = 0; nm < num; nm++)
    {
        GET_SCIA_ROE_INFO_MOD(eclipseMode, julianDay[nm],
                              orbit+nm, &saaFlag, orbitPhase+nm);
        //printf("orb %d = %d\n", nm, orbit[nm]);
    }
}
