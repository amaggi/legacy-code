/* includes */
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

/* defines */
#define CMT 1 /* catalogue types */ 
#define PDE 2
#define EHB 3
#define TXT 4
#define JDY 5

#define MOM 1 /* mag search type */
#define MAG 2
#define MS  3
#define MB  4

#define OUT_TXT 1 /* output format */
#define OUT_JDY 2 
#define OUT_REG 3 

#define SHORT 1 /* header format */
#define LONG 2 
#define NONE 3 

#define SET 1
#define UNSET 0

/* structures */
struct Event {
  struct tm origin_time;
  int msec; // just keeps the miliseconds - does no calculations on them
  float lat, lon, dep, mb, ms;
  float cmt_lat, cmt_lon, cmt_dep, mom_base;
  int strike1, dip1, rake1, strike2, dip2, rake2, mom_exp; 
  char region[80];
};

struct JDY_Record{
  char year[4], space1[2], jday[3], space2[2], month[2], space3[1],
       day[2], space4[1], hour[2], space5[1], min[2], space5a[1], sec[4],
       space6[2], lat[7], space7[2], lon[8], space8[2], depth[3],
       space9[2], mag[3], space10[1], name[13]; 
};

struct EHB_Record {
  char ahyp[1], isol[3], iseq[2], yr[2], space1[1], mon[2], space2[1]; 
  char day[2], space3[2], hr[2], space4[1], min[2];
  char sec[6], ad[1], lat[8], lon[8], dep[6], dmag[4];
  char imag[2], dummy1[62]; 
};

struct PDE_Record {
   char   blank1[1];
   char   Catalog_Source[5], blank2[1], Year[4], blank3[1], Month[2], Day[2], Hour[2], Min[2], Sec[5];
   char   Auth[2], Lat[7], Lon[8], Depth[3], dummy[2], DCD[1], pP_Phases[2], Stdev[4]; 
   char   mb[3], mbobs[2], Ms[3], ZH[1], Msobs[2], Mag1[4], Mag1_Scale[2];
   char   Mag1_Donor[5], Mag2[4], Mag2_Scale[2], Mag2_Donor[5], Reg[3];
   char   Sta[3], Io[1], Cult_effect[1], Isoseismal[1], FPS[1], Mom_Flag[1];
   char   ISC_Depth_Flag[1], IDE_Flag[1], Preferred_Flag[1], Flag[1];
   char   Phenomena_Codes[7], Radial_Dist[7];
};

struct CMT_Record{
  /*line 1*/
  char ident[8], space1[1], month[2], slash1[1], day[2], slash2[1], two_digit_year[2];
  char space2[1], hour[2], colon1[1], min[2], colon2[1], second[4], space3[1];
  char lat[6], space4[1], lon[7], space5[1], depth[5], mb[3], ms[3], region[24];
  char newline1[1];
  /* line 2*/
  char dummy1[44], cmt_lat[6], dummy2[6], cmt_lon[7], dummy3[6], cmt_dep[5], dummy4[5], newline2[1];
  /* line 3*/
  char dummy5[12], cmt_mom[2], dummy6[66], newline3[1];
  /* line 4*/
  char dummy7[45], cmt_base[4], space6[1], strike1[3], sapce2[1], dip1[2], space7[1], rake1[4], space8[1], strike2[3], space9[1], dip2[2], space10[1], rake2[4];
};

/* global variables */
int cat_type, cmt_hypo, mag_search, dep_search, circ_search, seg_search, rec_search, date_search, output_format, header_format;
float minmo, maxmo, minmag, maxmag, mindep, maxdep, minlat, minlon, maxlat, maxlon, cenlat, cenlon, rmin, rmax, bazmin, bazmax;
char cat_filename[FILENAME_MAX], basename[FILENAME_MAX];
struct tm start_time, end_time;
time_t mintime, maxtime;
FILE *fpin, *fpout;
struct Event event;


/* prototypes */

int parser (int argc, char *argv[]);
void print_help ();
void print_parsed_info();
void print_headers (int argc, char *argv[]);
int read_event();
void read_cmt_record();
void assign_cmt_event();
void read_pde_record();
void assign_pde_event();
void read_ehb_record();
void assign_ehb_event();
void read_jdy_record();
void assign_jdy_event();
void write_event();
int check_criteria();

float delta_angle(float a, float b);
float sum_angle(float a, float b);
extern void distaz_(float *elat, float *elon, float *slat, float *slon, float *azm, float *bzm, float *ddg, float *dkm);
extern int julday_(int *month, int *day, int *year);

void set_tm_jday(struct tm *tm);


int fixed_len_atoi(const char *string, int len);
double fixed_len_atof(const char *string, int len);
void copy_to_scratch_and_terminate(const char *string, int len);
