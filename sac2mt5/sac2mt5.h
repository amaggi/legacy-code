/* 
SAC2MT5
	header file

		* includes standard input and output macros and string 
		manipulation macros

		* defines constants

		* defines structures and turns them into type-definitions
			Compl
			Dat
			Npair
			Pz
			Sac
			Event

		* declares functions 
			resp2pz
			read_sac
			which_pz
			write_dsn
			which_instrument
			which_component
			read_txt
			write_header
			
Alessia Maggi	Completed March 1999
*/

#include<stdio.h>
#include<string.h>

#define MAXPOLES 50	/* maximum number of poles */
#define MAXZEROS 50	/* maximum number of zeros */
#define LNAMES   10	/* maximum length of names eg. station names */
#define FLEN     50	/* maximum length of filenames */  /* 40-> 50  Onur TAN,19tem02  */
#define MAXRESP  20	/* maximum number of responses per response file */
#define LINE    120	/* Maximum length of a line in a response file */
#define MAXPTS	5000	/* Maximum number of points in a seismogram */
#define MAXSTA	300	/* Maximum number of data files */
#define OK	0	
#define NOTOK	-1

typedef struct compl{		/* Structure for complex numbers */
        float re,im;
} Compl;

typedef struct dat{		/* Structure for data points */
	float x,y;
} Dat;

typedef struct name_pair{	/* Structure for pairs of data and response 
				filenames */
	char data[FLEN], resp[FLEN];
} Npair;

typedef struct pz{		/* Structure for polezero information */
        float A0, sensitivity;
        int n_poles, n_zeros;
        Compl poles[MAXPOLES];
        Compl zeros[MAXZEROS];
        int year, jday;
        char sta[LNAMES],netwk[LNAMES],channel[LNAMES];
        char fname[FLEN];
} Pz;

typedef struct sac{		/* Structure for Sac information */
	int npts, leven;
	float beg, ed, delta;
	int nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec;
	float origin;
	char staname[LNAMES], component[LNAMES], network[LNAMES];
	float stla, stlo, evla, evlo;
	Dat series[MAXPTS];
} Sac;

typedef struct event{		/* Structure for event information */
	int year, jday, month, day, hour, min, dep;
	float sec, lat, lon, mag;
} Event;


/* 	declaration of functions 	*/

int resp2pz(Pz poze[], FILE *fp);
int read_sac(Sac *data, FILE *fp);
int which_pz(Sac *data, Pz poze[], int n_pz);
int write_dsn(FILE *fp, Sac *data, Pz *poze);
int which_instrument(Sac *data);
int which_component(char *s, Sac *data);
int read_txt(FILE *fp, Event *quake);
int write_header(FILE *fp, Event *quake);

