/*
* $Id: evresp.h,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
* $Log: evresp.h,v $
* Revision 1.1.1.1  2002/07/12 11:15:19  maggi
*
*
* Revision 1.1  2002/05/23 10:28:52  maggi
* Initial revision
*
*
*/
/*==================================================================
 *                    evresp.h
 *=================================================================*/
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#define REVNUM "1.4"

#define PI (double)    3.141592654
#define TWOPI (double) 6.283185307
#define TRUE 1
#define FALSE 0
#define DELIM1 " \t\n"
#define STALEN 7
#define CHALEN 7
#define DATIMLEN 23
#define UNITSLEN 20


enum { UNDEF_UNITS, DIS, VEL, ACC };

/*  Filter types  */

enum { UNDEF_FILT, ANALOG, LAPLACE,
       FIR_SYM_1, FIR_SYM_2, FIR_ASYM, IIR_PZ
};

struct filter {
    int     lineno;
    int     filnum;
    int     type;
    int     nnum;
    int     nden;
    double  sint;
    double  ao;
    double  gain;
    double  *num;
    double  *den;
};

struct channel {
    char staname[STALEN];
    char chaname[CHALEN];
    char beg_t[DATIMLEN];
    char end_t[DATIMLEN];
    double gain;
    double srate;
    double sensit;
    double sensfreq;
    double calc_sensit;
    char inp_units[UNITSLEN];
    int units_code;
    char out_units[UNITSLEN];
    int seed_nstage;
    int nfilter;
    struct filter **filter;
    double del;
    double corr;
    double exp_del;
};


struct station {
    char name[STALEN];   
    char owner[40];
    char lat[20];     /* latitude,  decimal degress, + => N */
    char lon[20];     /* longitude, decimal degress, + => E */
    char elev[20];    /* elevation asl, meters              */
    char beg_t[DATIMLEN];
    char end_t[DATIMLEN];
    int nchannel;
    struct channel **channel;  /* array of pointers to channels */
};
