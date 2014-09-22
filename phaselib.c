#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <defs_nadc.h>
#include <defs_scia.h>
#include <proto_scia.h>

bool Use_Extern_Alloc = false;

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

void _GET_SCIA_ROE_ORBITPHASE(bool eclipseMode, bool saaFlag, int32_t num, int32_t orbit, float *orbitPhase, double *julianDay)
{
    int32_t nm;
    for (nm = 0; nm < num; nm++)
    {
        GET_SCIA_ROE_INFO( eclipseMode, julianDay[nm],
                           &orbit, &saaFlag, orbitPhase+nm );
    }
}

