/*
* $Id: evresp.c,v 1.1.1.1 2002/07/12 11:15:19 maggi Exp $
* $Log: evresp.c,v $
* Revision 1.1.1.1  2002/07/12 11:15:19  maggi
*
*
* Revision 1.1  2002/05/23 10:28:51  maggi
* Initial revision
*
*
*/
/*===================================================================
Name:      evresp Version 1.3
Purpose:
        Extract channel response parameters from either ASCII
        files produced by rdseed -s ("all-sta-cha" file) or
        rdseed -d ("sta-cha" files) and calculate the complex
        response.
Reference:
        SEED. Standard for the Exchange of Earthquake Data
        Reference Manual
        SEED Format Version 2.3
        February 1993
Author:    Jean-Francois Fels
Revisions:
    26 Jan 1994. Version 1.3
      - The name of the module is now evresp_(), callable by
        either a FORTRAN or a C program.
      - The 'output' array is now an array of single precision
        real numbers, declared in the calling program instead of 
        here.
    11 April 1994. Version ???
      - modified to handle a specified directory name in addition
        to reading from a specified file.
      - converted to unix-standard type command line inputs
		Thomas J. McSweeney:  tjm@iris.washington.edu

 *=================================================================*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "evresp.h"

#define NOERROR    0
#define NOSTATION -1
#define NOCHANNEL -2
#define NODATIME  -3
#define MAXENTRY 1

struct channel *Pchan[MAXENTRY];
char line[256];
int lineno = 0;
int err_line, err_lineno;
char *token[20];
char ifname[50];
char entry[MAXENTRY][50];
int nentry;
int pr_flag;
extern char *m_alloc();
int fil_count;                   /* unused in this code */



/*=================================================================
 *        evresp_()  C or FORTRAN callable module
 ================================================================*/
evresp_(sta,cha,datime,freq,output,nout,units,file,opt)
char    *sta;          /* station name */
char    *cha;          /* channel code */
char    *datime;       /* date: YYYY,DDD,HH:MM:SS */
float   *freq;         /* frequency */
float   *output;       /* n pairs real-imaginary */
int     *nout;         /* number of pairs */
char    *units;        /* units of response "acc"/"vel"/"dis"/null */
char    *file;
char    *opt;
{
    FILE *fp;
    struct channel *chan;
    int i, j, err;
    char *req_date;
    char req_args[100];
    int sta_found = FALSE;
    int cha_found = FALSE;
    int dat_found = FALSE;
    static char prev_args[100];
    static double out[MAXENTRY*2];
    struct stat buf;


    ucase(sta);
    ucase(cha);
    if (units != NULL) ucase(units);
    pr_flag = FALSE;
    if (opt != NULL && !strncmp(opt, "-v", 2)) {
        pr_flag = 1;
        if (opt[2] == '2') pr_flag = 2;
    }
    req_date = datime;

/* IF ARGS NOT NEW, CHANNEL ALREADY LOADED. GO COMPUTE RESPONSE */

    sprintf(req_args, "%s_%s_%s",sta,cha,req_date);
    if (!strncmp(req_args, prev_args, strlen(req_args)))
        goto channel_loaded;

/*==============    LOAD CHANNEL RESPONSE ================*/

    if (pr_flag) {
        fprintf(stderr,"<< IRIS SEED Channel Response ");
        fprintf(stderr,"Version 1 Release %s >>\n", REVNUM);
    }

/* STORE INPUT ARGS */

    sprintf(prev_args, "%s", req_args);

/* OPEN RESPONSE ASCII FILE */

    if ((err = open_file(file, sta, cha, ifname, &fp)))
        return err;

    chan = (struct channel*)m_alloc(0L,
            sizeof(struct channel), "chan");
    nentry = 0;
    Pchan[nentry] = chan;

/* LOOP THROUGH STATIONS-CHANNELS-DATE */

    if (file != NULL && strlen(file) != 0){
    	stat(file,&buf);
	if(!S_ISDIR(buf.st_mode))
		goto read_single_file;
    }

/*============= READ rdseed option -d FILES ================*/

    while (!(err = getline(fp, line, 255, 0, &lineno))) {
        if (!strncmp(line, "Station:", 8)) {

            sparse(line, token, DELIM1, 10);
            if (strncmp(token[1], sta, sizeof(sta))) continue;
            sprintf(chan->staname, "%s", sta);
            sta_found = TRUE;

            gettoken(fp,2,NULL,NULL);
            if (strncmp(token[1], cha, sizeof(cha))) continue;
            sprintf(chan->chaname, "%s", cha);
            cha_found = TRUE;

            get_date(fp,chan->beg_t,chan->end_t);
            if (cmp_date(chan->beg_t,chan->end_t, req_date)) {
                dat_found = TRUE;
                if (nentry >= MAXENTRY) break;
                sprintf(entry[nentry], "%s %s %s %s",
                                 sta,cha,chan->beg_t,chan->end_t);
                load_channel(fp, chan, stdout);
                Pchan[nentry++] = chan;
                chan = (struct channel*)m_alloc(0L,
                        sizeof(struct channel), "chan");
            }
        }
    }

    goto end_read;

read_single_file:

/*============= READ rdseed option -s FILE ================*/

    while (!(err = getline(fp, line, 255, 0, &lineno))) {

        if (!strncmp(line, "Station code:", 13)) {
            sparse(line, token, DELIM1, 10);
            sprintf(chan->staname, "%s", token[2]);
        }

        else if (!strncmp(line, "Channel:", 8)) {
            sparse(line, token, DELIM1, 10);
            sprintf(chan->chaname, "%s", token[1]);

            if (strncmp(sta,chan->staname,sizeof(chan->staname)))
                continue;
            sta_found = TRUE;

            if (!strncmp(cha,chan->chaname,sizeof(chan->chaname))) {
                cha_found = TRUE;
                while (!(err = getline(fp, line, 255, 0, &lineno)))
                    if (!strncmp(line, "Channel flags:", 14)) break;
                get_date(fp,chan->beg_t,chan->end_t);
                if (cmp_date(chan->beg_t,chan->end_t, req_date)) {
                    dat_found = TRUE;
                    load_channel(fp, chan, stdout);
                    Pchan[nentry++] = chan;
                    break;
                }
            }
        }
    }

end_read:

    if (!sta_found) return NOSTATION;
    if (!cha_found) return NOCHANNEL;
    if (!dat_found) return NODATIME;

    if (nentry > 1) {
        printf("===========================================\n");
        printf("WARNING: multiple entries.\n");
        if (nentry >= MAXENTRY) {
            nentry = MAXENTRY;
            printf("We'll keep only the %d following:\n",nentry);
        }
        for (i=0; i<nentry; i++) printf("%s\n",entry[i]);
        printf("===========================================\n");
    }

channel_loaded:

/*==============    CALCULATE  RESPONSE ================*/

    for (i=0,j=0; i<nentry; i++,j+=2) {
        calc_resp(Pchan[i], *freq, &out[i*2], units);
        output[j+0] = out[j+0];
        output[j+1] = out[j+1];
    }

    *nout  = nentry;
    return NOERROR;
}


/*=================================================================
 *                      open_file()
 ================================================================*/
open_file(req, sta, cha, ifname, fp)
char *req, *ifname;
FILE **fp;
{
char *basedir;
struct stat buf;

    if (req != NULL && strlen(req) != 0) {
        stat(req,&buf);
        if(S_ISDIR(buf.st_mode)){
            sprintf(ifname, "%s/RESP.%s.%s",req,sta,cha);
	}
        else{
            sprintf(ifname, "%s", req);
	}
        if ((*fp = fopen(ifname, "r")) == NULL) {
            fprintf(stderr,"evresp: can't open '%s'\n", ifname);
            return -4;
        }
        return 0;
    }
    sprintf(ifname, "RESP.%s.%s", sta,cha);
    if ((*fp = fopen(ifname, "r")) == NULL) {
        if ((basedir = (char *) getenv("SEEDRESP")) == NULL) {
            fprintf(stderr,"can't get env SEEDRESP\n");
            return -4;
        }
        sprintf(ifname, "%s/RESP.%s.%s", basedir,sta,cha);
        if ((*fp = fopen(ifname, "r")) == NULL) {
            fprintf(stderr,"evresp: can't open '%s'\n", ifname);
            return -4;
        }
    }
    return 0;
}

/* %W% %G% */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#define DELIM1 " \t\n"
#define COEFSIZE 32
#define FIR_NORM_TOL 0.02

extern char line[];
extern int lineno;
extern int err_line, err_lineno;
extern char *token[];
extern int fil_count;
extern char ifname[];
extern int pr_flag;
extern char *m_alloc();


/*=================================================================
 *                   load_channel() 
 ================================================================*/
load_channel(fp, chan, ofp)
FILE *fp, *ofp;
struct channel *chan;
{
int err, i, j;
int flt_count;
struct filter *filt1, *prev_filt;
int stage_num, prev_stage_num;
int filt_continuation;
int no_filter;
int new_filter;
char resp_type[100], tf_type[100];
int lineno_filt;
double out[2];
char units[UNITSLEN];

    flt_count = 0;
    filt1 = NULL;
    chan->del = 0.0;
    chan->corr = 0.0;
    chan->gain = 1.0;
    chan->seed_nstage = 0;
    stage_num = prev_stage_num = -1;
    filt_continuation = FALSE;
    no_filter = FALSE;

    while (!(err = getline(fp, line, 255, 0, &lineno))) {

        if (!strncmp(line, "Response type:", 14)) {
            sparse(line, token, DELIM1, 10);
            sprintf(resp_type, "%s", token[2]);
            lineno_filt = lineno;

            if (!strncmp(resp_type,"Poles",5) ||
                !strncmp(resp_type,"Coefficients",12)) {

        /* GET TF TYPE, STAGE NUMBER, UNITS */

                prev_stage_num = stage_num;

                gettoken(fp, 3, "Transfer function type:",tf_type);
                gettoken(fp, 3, "Stage sequence",NULL);
                stage_num = atoi(token[3]);
                chan->seed_nstage = stage_num;

                if (prev_stage_num == -1 && stage_num == 0) {
                    printf("WARNING: First stage in sequence is ");
                    printf("null; %s %s %s\n", 
                        chan->staname, chan->chaname, chan->beg_t);
                }

                read_units(fp,chan,resp_type,tf_type,stage_num);

        /* CREATE A NEW FILTER */

                mem("filter", &chan->filter, flt_count);
                filt1 = chan->filter[flt_count++];
                chan->nfilter = flt_count;

                if (stage_num == prev_stage_num)
                    filt_continuation = TRUE;

        /* READ PZ FILTER */
                if (!strncmp(resp_type, "Poles", 5))
                    new_filter = read_pzfilt(fp, filt1, tf_type);

        /* READ COEF FILTER */
                else if (!strncmp(resp_type, "Coefficients", 12))
                    new_filter = read_coefilt(fp,filt1, tf_type,
                                     chan,filt_continuation);

                filt_continuation = FALSE;

        /* IF NO TRUE FILTER FOUND, DEALLOCATE MEMORY */

                if (new_filter == FALSE) {
                    mem("filter", &chan->filter, --flt_count);
                    filt1 = chan->filter[flt_count-1];
                    chan->nfilter = flt_count;
                }

            } /* end if resp_type == Poles || Coefficients */

        /* READ DECI_FAC BLOCKETTE */

            else if (!strncmp(resp_type, "Decimation", 10)) {
                prev_stage_num = stage_num;
                gettoken(fp,3,"Stage sequence",NULL);
                stage_num = atoi(token[3]);
                chan->seed_nstage = stage_num;

                if (filt1 == NULL || stage_num != prev_stage_num)
                    no_filter = TRUE;

                read_decim(fp, filt1, chan, no_filter);

                no_filter = FALSE;

            } /* end if resp_type == Decimation */

        } /* end if line == "Response type" */

    /* READ GAIN BLK */

        else if (!strncmp(line, "Stage sequence number:", 22)) {
            sparse(line, token, DELIM1, 10);
            stage_num = atoi(token[3]);
            if (stage_num != 0) chan->seed_nstage = stage_num;

            read_gain(fp, chan, stage_num);

    /* END OF CHANNEL */

            if (stage_num == 0) break;

        } /* end if line == "Stage sequence number" */
    }

/* MAKE SYMETRICAL FILTERS IF POSSIBLE */

    for (j=0; j<chan->nfilter; j++) {
        if (chan->filter[j]->type == FIR_ASYM)
            check_sym(chan->filter[j], chan);
    }

    fix_channel(chan);

/* COMPUTE SENSITIVITY AT SENS_FREQ FOR CHECKING */

    if      (chan->units_code == DIS) sprintf(units, "DIS");
    else if (chan->units_code == VEL) sprintf(units, "VEL");
    else if (chan->units_code == ACC) sprintf(units, "ACC");
    else                              sprintf(units, "");
    calc_resp(chan, chan->sensfreq, out, units);
    chan->calc_sensit = sqrt(out[0]*out[0]+out[1]*out[1]);

    if (pr_flag) {
        print_channel(ofp, chan);
        fprintf(ofp, "\n");
    }
}



/*=================================================================
 *                      read_units()
 *================================================================*/
read_units(fp, chan, resp_type, tf_type, stage_num)
FILE *fp;
struct channel *chan;
char *resp_type, *tf_type;
int stage_num;
{
char inp_units[UNITSLEN];

    gettoken(fp,5,"Response in",inp_units);
    gettoken(fp,5,"Response out",chan->out_units);

   /* IF SENSOR, SAVE INPUT UNITS */

    if (!strncmp(resp_type,"Poles", 5) &&
        (!strncmp(tf_type,  "Analog",6) ||
         !strncmp(tf_type,  "Laplace",7))) {
        if (stage_num == 1)
            write_units(inp_units, chan->inp_units, &(chan->units_code));
    }
}


/*=================================================================
 *                      read_gain()
 *================================================================*/
read_gain(fp, channel, stage_num)
struct channel *channel;
FILE *fp;
int stage_num;
{
static int first = TRUE;

    if (first) {
        first = FALSE;
        channel->gain = 1.0;
    }
    if (stage_num != 0) {
        gettoken(fp,1,"Gain",NULL);
        channel->gain *= atof(token[1]);
        skiplines(fp, 1);
    }
    else {
        gettoken(fp,1,"Sensitivity",NULL);
        channel->sensit = atof(token[1]);
        gettoken(fp,4,"Frequency",NULL);
        channel->sensfreq = atof(token[3]);
    }
    skiplines(fp, 1);
    getline(fp, line, 255, 0, &lineno);
    if (strlen(line) != 0) {
       fprintf(stderr,"ERROR (read_gain): unexpected end of data\n");
       fprintf(stderr,"%s line %d\n", ifname, lineno);
       fprintf(stderr,"\tExecution terminating.\n");
       exit(1);
    }
}


/*=================================================================
 *                       read_decim()
 *================================================================*/
read_decim(fp, filter, channel, no_filter)
FILE *fp;
struct filter *filter;
struct channel *channel;
int no_filter;
{
double srate;
int decim;
char del_a[20];
char corr_a[20];
static first = TRUE;
struct filter *prev_filt;
int fnum, nc;

    if (first == TRUE) {
        channel->del  = 0.0;
        channel->corr = 0.0;
        channel->exp_del = 0.0;
        channel->srate = 0.0;
        first = FALSE;
    }

/* GET INP_SRATE and DECIFAC */

    gettoken(fp,4,"Input sample rate",NULL);
    srate = atof(token[3]);
    gettoken(fp,2,"Decimation factor",NULL);
    decim = atoi(token[2]);
    channel->srate = srate / (double)decim;

/* WRITE INP_SAMPLE INTERVAL INTO FILTER */

    if (no_filter == FALSE) {
        if (srate != 0.0) filter->sint = 1.0/srate;
        else              filter->sint = 0.0;
    }

    skiplines(fp, 1);

/* GET FILTER DELAY AND CORRECTION */

    gettoken(fp,3,"Estimated delay (seconds)",del_a);
    gettoken(fp,3,"Correction applied (seconds)",corr_a);
    channel->del  += atof(del_a);
    channel->corr += atof(corr_a);

    getline(fp, line, 255, 0, &lineno);
    if (strlen(line) != 0) {
      fprintf(stderr,"ERROR (read_decim): unexpected end of data\n");
      fprintf(stderr,"%s line %d\n", ifname, lineno);
      fprintf(stderr,"\tExecution terminating.\n");
      exit(1);
    }
}



/*=================================================================
 *                      read_pzfilt()
 *================================================================*/
read_pzfilt(fp, f, type)
FILE *fp;
struct filter *f;
char *type;
{
int i, j, nnum, nden;
char pattern[6];


    if      (!strncmp(type, "Analog", 6))  f->type = ANALOG;
    else if (!strncmp(type, "Laplace", 7)) f->type = LAPLACE;
    else if (!strncmp(type, "Digital", 7)) f->type = IIR_PZ;
    else {
	fprintf(stderr, "ERROR (read_pzfilt): filter type ");
	fprintf(stderr, "'%s' not supported\n", type);
	fprintf(stderr,"%s line %d\n",ifname, lineno);
	fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    f->sint = 0.0;
    f->gain = 1.0;

    gettoken(fp,3,"AO norm",NULL);      f->ao = atof(token[3]);
    skiplines(fp, 1);
    gettoken(fp,3,"Number of z",NULL);  f->nnum = atoi(token[3]);
    gettoken(fp,3,"Number of p",NULL);  f->nden = atoi(token[3]);

    nnum = f->nnum;
    nden = f->nden;
    if (nnum == 0 && nden == 0)
        return FALSE;

    if (nnum) f->num = (double *)m_alloc(0L,
                        nnum*sizeof(double)*2, "f->num");
    if (nden) f->den = (double *)m_alloc(0L,
                        nden*sizeof(double)*2, "f->den");

    gettoken(fp,0,"Complex z",NULL);
    skiplines(fp, 1);
    for (i=0,j=0; i<nnum; i++,j+=2) {
        sprintf(pattern, "%3d", i);
        gettoken(fp,5,pattern,NULL);
        f->num[j]   = atof(token[1]);
        f->num[j+1] = atof(token[2]);
    }

    gettoken(fp,0,"Complex p",NULL);
    skiplines(fp, 1);
    for (i=0,j=0; i<nden; i++,j+=2) {
        sprintf(pattern, "%3d", i);
        gettoken(fp,5,pattern,NULL);
        f->den[j]   = atof(token[1]);
        f->den[j+1] = atof(token[2]);
    }

    getline(fp, line, 255, 0, &lineno);
    if (strlen(line) != 0) {
       fprintf(stderr,"ERROR (read_pzf): unexpected end of data\n");
       fprintf(stderr,"%s line %d\n", ifname, lineno);
       fprintf(stderr,"\tExecution terminating.\n");
       exit(1);
    }

    return TRUE;
}


/*=================================================================
 *                      read_coefilt()
 *================================================================*/
read_coefilt(fp, f, type, chan, filt_continuation)
FILE *fp;
struct filter *f;
char *type;
struct channel *chan;
int filt_continuation;
{
int i, j, nnum, nden;
char pattern[6];
struct filter *prev_f;

    if (!strncmp(type,"D", 1)) f->type = FIR_ASYM;
    else {
        fprintf(stderr, "ERROR (read_coefilt): filter type ");
        fprintf(stderr, "'%s' not supported\n", type);
        fprintf(stderr,"%s line %d\n",ifname, lineno);
	fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    f->sint = 0.0;
    f->ao   = 1.0;
    f->gain = 1.0;

    gettoken(fp,3,"Number of num",NULL);
    f->nnum = atoi(token[3]);
    if ((nnum = f->nnum) == 0) return FALSE;

    gettoken(fp,3,"Number of den",NULL);
    f->nden = atoi(token[3]);
    if ((nden = f->nden) != 0) {
        fprintf(stderr,"ERROR (read_coefilt): n of denominators ");
        fprintf(stderr,"not null;  not supported filter type\n");
        fprintf(stderr,"%s line %d\n",ifname, lineno);
        fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }

    f->num = (double *)m_alloc(0L,nnum*sizeof(double), "f->num");
    f->den = NULL;

    if (nnum != 0) {
        gettoken(fp,0,"Numerator",NULL);
        skiplines(fp, 1);
    }
    for (i=0; i<nnum; i++) {
        sprintf(pattern, "%3d", i);
        gettoken(fp,3,pattern,NULL);
        f->num[i] = atof(token[1]);
    }

    getline(fp, line, 255, 0, &lineno);
    if (strlen(line) != 0) {
       fprintf(stderr,"ERROR (read_coef): unexpected end of data\n");
       fprintf(stderr,"file %s line %d\n", ifname, lineno);
       fprintf(stderr,"\tExecution terminating.\n");
       exit(1);
    }

    if (filt_continuation == TRUE) {
        prev_f = chan->filter[chan->nfilter -2];
        if (prev_f->type != FIR_ASYM) {
            fprintf(stderr,"ERROR (read_coefilt): found multiple ");
            fprintf(stderr,"filters of different types \n");
            fprintf(stderr,"in %s %s stage %d\n",
                chan->staname, chan->chaname, chan->seed_nstage);
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
        i = prev_f->nnum;
        prev_f->nnum += f->nnum;
        prev_f->num = (double *)m_alloc(prev_f->num,
                       prev_f->nnum*sizeof(double), "prev_f->num");
        for (j=0; i<prev_f->nnum; i++, j++) {
            prev_f->num[i] = f->num[j];
        }
        return FALSE;
    }

    return TRUE;
}



/*=================================================================
 *                       get_date()
 *================================================================*/
get_date(fp, beg, end)
FILE *fp;
char *beg, *end;
{
char beg_time[12], end_time[12];
int beg_yd, end_yd;
int beg_t, end_t;

    if (beg == NULL || end == NULL) {
        fprintf(stderr,"ERROR (get_date): null pointer in args\n");
	fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }

    sprintf(beg, "XXXX,XXX,00:00:00");
    gettoken(fp, 2, "Start date:",NULL);
    strncpy(beg, token[2], strlen(token[2])); beg[17] = 0;

    sprintf(end, "XXXX,XXX,00:00:00");
    gettoken(fp, 2, "End date:",  NULL);
    strncpy(end, token[2], strlen(token[2])); end[17] = 0;
    if (!strncmp(token[2], "No",2)) sprintf(end,"2041,001,00:00:00");
    if (!strncmp(token[2], "(n",2)) sprintf(end,"2041,001,00:00:00");
}



/*===================================================================
 *  Read a single line from the given file, stripping out comments
 *  and blank lines.
 *  The processed line will be a NULL terminated string and without
 *  the trailing newline.
 *  Return values: 0 => success
 *                 1 => EOF
 *                 2 => read or other error
 *=================================================================*/
int getline(fp, buffer, buflen, comment, lineno)
FILE *fp;     /* input stream                   */
char *buffer; /* buffer to hold line            */
int  buflen;  /* length of buffer               */
char comment; /* comment character, not handled */
int  *lineno; /* line number of line in buffer  */
{
int i;

    clearerr(fp);
    buffer[0] = 0;

/* READ NEXT LINE IN THE FILE */

    if (fgets(buffer, buflen-1, fp) == NULL) {
        buffer[0] = 0;
        if (feof(fp))        return 1;
        else {
            err_line = TRUE;
            err_lineno = ++*lineno;
            return 2;
        }
    }
    ++*lineno;

/* REMOVE TRAILING BLANKS */

    i = strlen(buffer) - 1;
    while (i >= 0 && (buffer[i] == ' '|| buffer[i] == '\n')) --i;
    buffer[++i] = 0;
        
    return 0;
}



/*=================================================================
 *                 Get tokens
 *=================================================================*/
gettoken(fp, toknum, keyword, dest)
FILE *fp;
char *keyword, *dest;
int toknum;
{
int ntok;

    getline(fp, line, 255, 0, &lineno);
    if (keyword != NULL && strncmp(line,keyword,strlen(keyword))) {
        fprintf(stderr,"ERROR (gettoken): pattern %s ", keyword);
        fprintf(stderr,"not found lineno %d\n", lineno);
	fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    ntok = sparse(line, token, DELIM1, 10);
    if (ntok < toknum) {
        fprintf(stderr,"ERROR (gettoken): found only %d ", ntok);
        fprintf(stderr,"tokens %s lineno %d\n", ifname, lineno);
	fprintf(stderr,"\tExecution terminating.\n");
	exit(1);
    }
    if (dest != NULL) sprintf(dest, "%s", token[toknum]);

    return ntok;
}



/*=================================================================
 *                  Skip lines
 *=================================================================*/
skiplines(fp, n)
FILE *fp;
int n;
{
    int i;
    for (i=0; i<n ; i++) {
        if (fgets(line, 255, fp) == NULL) {
            fprintf(stderr, "ERROR. (skiplines): ");
	    fprintf(stderr, " %s line %d\n",ifname, lineno);
	    fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
    }
    lineno += n;
}



/*==================================================================
 *       Allocate, reallocate memory for variable length lists
 *=================================================================*/
mem(type, list, n)
char *type;
char ***list;
int n;
{
struct channel **chalist;
struct filter  **fltlist;

    if      (!strncmp(type, "channel", 7)) {
        chalist = (struct channel **) *list;
        chalist = (struct channel **) m_alloc(chalist,
                   sizeof(struct channel *) * (n+1), "channel_p");
        chalist[n] =(struct channel *)m_alloc(0L,
                   sizeof(struct channel), "channel");
        *list = (char **) chalist;
    }

    else if (!strncmp(type, "filter",  6)) {
        fltlist = (struct filter  **) *list;
        fltlist = (struct filter  **) m_alloc(fltlist,
                   sizeof(struct filter  *) * (n+1), "filter_pt");
        fltlist[n] =(struct filter  *)m_alloc(0L,
                   sizeof(struct filter ), "filter");
        *list = (char **) fltlist;
    }

    else {
        if (type && strlen(type)) fprintf(stderr,
	    "ERROR (mem): type '%s' not supported\n", type);
        else fprintf(stderr,
	    "ERROR (mem): arg 'type' missing\n");
        fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
}

/* %W% %G% */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>


#define FIR_NORM_TOL 0.02

extern int   pr_flag;
extern char *m_alloc();


/*=================================================================
 *                       print_channel()
 *=================================================================*/
print_channel(fp, ch)
FILE *fp;
struct channel *ch;
{
int j;

    fprintf(fp, "%s %s ", ch->staname, ch->chaname);
    fprintf(fp, "%s %s   Seed units %s/%s\n",
        ch->beg_t, ch->end_t, ch->inp_units, ch->out_units);
    fprintf(fp, "gain=%.5E sens=%.5E (calc=%.5E) @ %.5E Hz\n",
        ch->gain, ch->sensit, ch->calc_sensit, ch->sensfreq);
    fprintf(fp, "del=%.5E  corr=%.5E  exp=%.5E  srate=%.3fsps\n",
        ch->del, ch->corr, ch->exp_del, ch->srate);
    for (j=0; j<ch->nfilter; j++) {
        print_filter(fp, ch->filter[j], j);
    }
    fflush(fp);
}

/*=================================================================
 *                      print_filter()
 ================================================================*/
print_filter(fp, f, n)
FILE *fp;
struct filter *f;
int n;
{
int type;
char atype[12];

    type = f->type;
    if      (type == ANALOG)    sprintf(atype, "ANALOG   ");
    else if (type == LAPLACE)   sprintf(atype, "LAPLACE  ");
    else if (type == FIR_SYM_1) sprintf(atype, "FIR_SYM_1");
    else if (type == FIR_SYM_2) sprintf(atype, "FIR_SYM_2");
    else if (type == FIR_ASYM)  sprintf(atype, "FIR_ASYM ");
    else if (type == IIR_PZ)    sprintf(atype, "IIR_PZ   ");
    else                        sprintf(atype, ".........");

    if (n >= 0) fprintf(fp, "Filter %d ", n+1);
    fprintf(fp, "%s %E %E %E %3d %3d\n",
        atype, f->sint, f->ao, f->gain, f->nnum, f->nden);
}


/*=================================================================
 *                       fix_channel()
 *================================================================*/
fix_channel(ch)
struct channel *ch;
{
extern int nentry;
struct filter *f;
double del, corr;
int i, j, nc;

    del = 0.0; corr = 0.0;
    for (i=0; i<ch->nfilter; i++) {
        f = ch->filter[i];
        if (f->type == ANALOG || f->type == LAPLACE);
        else if (f->sint == 0.0) {
            printf("WARNING: FIR discarded: null sint; ");
            printf("%s %s\n", ch->staname, ch->chaname);
            for (j=i; j<ch->nfilter-1; j++)
                ch->filter[j] = ch->filter[j+1];
            --(ch->nfilter);
            --(ch->seed_nstage);
        }
    }

/* CALCULATE FILTER EXPECTED DELAY */

    ch->exp_del = 0;
    for (i=0; i<ch->nfilter; i++) {
        f = ch->filter[i];
        if      (f->type == FIR_SYM_1) nc = f->nnum*2 + 1;
        else if (f->type == FIR_SYM_2) nc = f->nnum*2;
        else if (f->type == FIR_ASYM)  nc = f->nnum;
        else nc = 0;
        if(nc) ch->exp_del += ((nc-1)/2.0) * f->sint;
    }

    if (ch->seed_nstage == 1)
        ch->gain *= ch->sensit;

    if (ch->gain == 0.0) {
        printf("WARNING: %s.%s %s : null gain, set to 1.\n",
            ch->staname, ch->chaname, ch->beg_t);
        ch->gain = 1.0;
    }
}



/*=================================================================
 *                       check_sym()
 *================================================================*/
check_sym(f, chan)
struct filter *f;
struct channel *chan;
{
int nc, n0, k;
double sum = 0.0;

    nc = f->nnum;

/* CHECK IF IF FILTER IS NORMALIZED TO 1 AT FREQ 0 */

    for (k=0; k<nc; k++) sum += f->num[k];
    if (sum < (1.0-FIR_NORM_TOL) || sum > (1.0+FIR_NORM_TOL)) {
        printf("WARNING: FIR normalized: sum[coef]=%E; ", sum);
        printf("%s %s\n", chan->staname, chan->chaname);
        for (k=0; k<nc; k++) f->num[k] /= sum;
    }

    if (f->type != FIR_ASYM) return;

/* CHECK IF FILTER IS SYMETRICAL WITH EVEN NUM OF WEIGHTS */

    if ((nc%2) == 0) {
        n0 = nc / 2;
        for (k=0; k < n0; k++) {
            if (f->num[n0+k] != f->num[n0-k-1]) return;
        }
        f->type = FIR_SYM_2;
        f->nnum = n0;
        for (k=0; k < n0; k++) f->num[k] = f->num[n0+k];
    }

/* CHECK IF FILTER IS SYMETRICAL WITH ODD NUM OF WEIGHTS */

    else {
        n0 = (nc - 1) / 2;
        for (k=1; k<nc-n0; k++) {
            if (f->num[n0+k] != f->num[n0-k]) return;
        }
        f->type = FIR_SYM_1;
        f->nnum =  nc-n0;
        for (k=0; k<nc-n0; k++) f->num[k] = f->num[n0+k];
    }
}


/*=================================================================
 *                       write_units()
 *================================================================*/
write_units(src_u, dest_u, dest_ucode)
char *src_u, *dest_u;
int *dest_ucode;
{
    sprintf(dest_u, "%s", src_u);
    if      (strlen(src_u)==6 && !strncmp(src_u, "M/S**2", 6))
        *dest_ucode = ACC;
    else if (strlen(src_u)==3 && !strncmp(src_u, "M/S",    3))
        *dest_ucode = VEL;
    else if (strlen(src_u)==1 && !strncmp(src_u, "M",      1))
        *dest_ucode = DIS;
}
/*=================================================================*/



/*=================================================================
 *                       cmp_date()
 *================================================================*/
cmp_date(beg, end, req)
char *beg, *end, *req;
{
char beg_time[12], end_time[12], req_time[12];
int beg_yd, end_yd, req_yd;
int beg_t, end_t, req_t;
int l;

    if (beg == NULL || end == NULL || req == NULL) {
        fprintf(stderr,"ERROR (cmp_date): null pointer in args\n");
        fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    if (!strlen(beg) || !strlen(end) || !strlen(req)) {
        fprintf(stderr,"ERROR (cmp_date): null date\n");
        fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    if (!strncmp(req, "anytime", 7)) return TRUE;

    beg_yd = atoi(&beg[0])*1000 + atoi(&beg[5]);
    end_yd = atoi(&end[0])*1000 + atoi(&end[5]);
    req_yd = atoi(&req[0])*1000 + atoi(&req[5]);

    if ((l=strlen(beg))<8 || beg_yd<1900001 || beg_yd>2099365) {
        fprintf(stderr,"ERROR (cmp_date): ");
        fprintf(stderr,"invalid beg date: %s\n", beg);
        fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    if ((l=strlen(end))<8 || end_yd<1900001 || end_yd>2099365) {
        fprintf(stderr,"ERROR (cmp_date): ");
        fprintf(stderr,"invalid end date: %s\n", end);
        fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
/*
    printf("   %d %d %d \n", beg_yd,end_yd,req_yd);
*/
    if (req_yd >= beg_yd && req_yd <= end_yd) {
        sprintf(beg_time, "00:00:00");
        sprintf(end_time, "00:00:00");
        sprintf(req_time, "00:00:00");

        strncpy(beg_time, &beg[9], strlen(&beg[9]));
        strncpy(end_time, &end[9], strlen(&end[9]));
        strncpy(req_time, &req[9], strlen(&req[9]));

        beg_t = atoi(&beg_time[0])*10000 +
                atoi(&beg_time[3])*100 +
                atoi(&beg_time[6]);
        end_t = atoi(&end_time[0])*10000 +
                atoi(&end_time[3])*100 +
                atoi(&end_time[6]);
        req_t = atoi(&req_time[0])*10000 +
                atoi(&req_time[3])*100 +
                atoi(&req_time[6]);
/*
ANMO EHZ 1990,012,00:00:00 1990,253,14:04:00
ANMO EHZ 1990,253,14:04:00 1991,025,00:00:00
    printf("   %d %d %d \n", beg_t,end_t,req_t);
*/
        if (req_yd == beg_yd) {
            if (beg_t != 0) {
                if (req_t == 0) {
                    fprintf(stderr,"ERROR: Found start time: %s. ",
                        beg);
                    fprintf(stderr,"Please specify time of day.\n");
                    exit(1);
                }
                if (req_t >= beg_t) return TRUE;
                if (req_t < beg_t) return FALSE;
            }
        }
        else if (req_yd == end_yd) {
            if (end_t != 0) {
                if (req_t == 0) {
                    fprintf(stderr,"ERROR: Found end   time: %s. ",
                        end);
                    fprintf(stderr,"Please specify time of day.\n");
                    exit(1);
                }
                if (req_t <= end_t) return TRUE;
                if (req_t > end_t) return FALSE;
            }
        }
        return TRUE;
    }
    return FALSE;
}



/*=================================================================
 *                   Make Upper Cases
 *=================================================================*/
ucase(c)
char *c;
{
int i, l;
    for (i=0;i<(l=strlen(c));i++)
        if (islower(c[i])) c[i] = toupper(c[i]);
}


/*===================================================================
 *                 Parse the given string.
 *=================================================================*/
int sparse(input, argv, delimiters, max_tokens)
char *input;
char *argv[];
char *delimiters;
int  max_tokens;
{
extern char *strtok();
int i = 0;

    if (max_tokens < 1) {
        fprintf(stderr,"ERROR (sparse): illegal 'max_tokens'\n");
	fprintf(stderr,"\tExecution terminating.\n");
        exit(1);
    }
    i = 0;
    if ((argv[i] = strtok(input, delimiters)) == NULL) return 0;
    for (i = 1; i < max_tokens; i++) {
        if ((argv[i] = strtok(NULL, delimiters)) == NULL) return i;
    }
    return i;
}


/*==================================================================
 *                  Allocate, reallocate memory
 *=================================================================*/
char *m_alloc(ptr, nbytes, caller)
char *ptr;
int nbytes;
char *caller;
{
char *p;
int dbug = FALSE;
FILE *fp = NULL;

    if (dbug) fp = stdout;

    if (ptr == NULL && nbytes == NULL) {
        return NULL;
    }
    if (nbytes < 0) {
        if (ptr == NULL) {
            fprintf(stderr,"m_alloc: cannot free a null pointer");
            if (caller && strlen(caller))
                 fprintf(stderr," ('%s')\n", caller);
            else fprintf(stderr,"\n");
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
        if (free(ptr) == NULL) {
            fprintf(stderr,"m_alloc: free failed");
            if (caller && strlen(caller))
                 fprintf(stderr," in '%s'\n", caller);
            else fprintf(stderr,"\n");
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
        return NULL;
    }

    if (ptr == NULL) {
        if ((p = (char *) malloc(nbytes)) == NULL) {
            fprintf(stderr,"m_alloc: malloc failed");
            if (caller && strlen(caller))
                 fprintf(stderr," in '%s'\n", caller);
            else fprintf(stderr,"\n");
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
        memset(p, 0, nbytes);
        if (dbug && caller && strlen(caller))
            fprintf(fp,"m_alloc: alloc   '%s' %d\n", caller, nbytes);
    }
    else {
        if((p = (char *) realloc(ptr, nbytes)) == NULL) {
            fprintf(stderr,"m_alloc: realloc failed");
            if (caller && strlen(caller))
                 fprintf(stderr," in '%s'\n", caller);
            else fprintf(stderr,"\n");
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
        if (dbug && caller && strlen(caller))
            fprintf(fp,"m_alloc: realloc '%s' %d\n", caller, nbytes);
    }
/*    if (dbug && fp) fflush(fp); */
    fflush(stdout);
    return p;
}



/*=================================================================
*                   Calculate response
*=================================================================*/
calc_resp(chan, f, output, out_units)
struct channel *chan;
double          f;
double          *output;
char            *out_units;
{
struct filter *fil;
int j;
double w, wsint, ao;
double *num, *den;
int nnum, nden;
double of[2];
double rpart, ipart;
double mod;

    w = TWOPI * f;
    rpart = 1.0; ipart =  0.0;

/* LOOP THRU FILTERS */

    for (j=0; j<chan->nfilter; j++) {

        fil = chan->filter[j];
        if (fil == NULL) {
            fprintf(stderr,"ERROR (calc_resp): null filter\n");
	    fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
        nnum = fil->nnum; num = fil->num;
        nden = fil->nden; den = fil->den;
        if (num == NULL && den == NULL) {
            fprintf(stderr,"ERROR (calc_resp): null pointer to ");
            fprintf(stderr,"filter\n");
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }

        ao = fil->ao;
        wsint = w * fil->sint;

        if (fil->type == ANALOG || fil->type == LAPLACE)
            analog_trans(num, nnum, den, nden, ao, f, of, fil->type);

        else if (fil->type == IIR_PZ)
            iir_pz_trans(num, nnum, den, nden, ao, wsint, of);

        else if (fil->type == FIR_SYM_1 || fil->type == FIR_SYM_2)
            fir_sym_trans(num, nnum, ao, wsint, of, fil->type);

        else if (fil->type == FIR_ASYM)
            fir_asym_trans(num, nnum, ao, w, fil->sint, of);

        else {
           fprintf(stderr,"ERROR (calc_resp): filter type ");
           fprintf(stderr,"%d not supported\n", fil->type);
           fprintf(stderr,"\tExecution terminating.\n");
           exit(1);
        }
        zmul(of[0], of[1], rpart, ipart, &rpart, &ipart);
        if (pr_flag == 2)
            printf("%f R=%.10e I=%.10e\n", f, of[0], of[1]);
    }

/*  WRITE REAL PART IN output[0], IMAG PART in output[1] */

    output[0] = rpart * chan->gain;
    output[1] = ipart * chan->gain;

    convert_to_units(chan->units_code, out_units, output, w);

    if (pr_flag == 2) {
        mod = sqrt(output[0]*output[0] + output[1]*output[1]);
        printf("%f        R=%.10e I=%.10e Mod=%.7e\n\n",
                 f,       output[0], output[1], mod);
    }
}


/*==================================================================
 * Convert response to velocity first, then to specified units
 *=================================================================*/
convert_to_units(inp, out_units, data, w)
int inp;
char *out_units;
double *data;
double w;
{
int out, l;

    if (out_units != NULL && (l=strlen(out_units)) > 0) {
        if      (!strncmp(out_units, "DIS", 3)) out = DIS;
        else if (!strncmp(out_units, "VEL", 3)) out = VEL;
        else if (!strncmp(out_units, "ACC", 3)) out = ACC;
        else {
            fprintf(stderr,"ERROR (convert_to_units): bad output ");
            fprintf(stderr,"units\n");
            fprintf(stderr,"\tExecution terminating.\n");
            exit(1);
        }
    }
    else out = VEL;

    if (inp == DIS) {
        if (out == DIS) return;
        if (w != 0.0)
            zmul(data[0],data[1], 0.0,-1.0/w, &data[0],&data[1]);
        else data[0] = data[1] = 0.0;
    }
    else if (inp == ACC) {
        if (out == ACC) return;
        zmul(data[0],data[1], 0.0,w, &data[0],&data[1]);
    }

    if (out == DIS) {
        zmul(data[0],data[1], 0.0,w, &data[0],&data[1]);
    }
    else if (out == ACC) {
        if (w != 0.0)
            zmul(data[0],data[1], 0.0,-1.0/w, &data[0],&data[1]);
        else data[0] = data[1] = 0.0;
    }

}


/*==================================================================
 *                Response of analog filter
 *=================================================================*/
analog_trans(ze, nz, po, np, h0, f, out, type)
double *ze, *po, h0, f, *out;
int nz, np;
int type;
{
    int i;
    double mod = 1.0, pha = 0.0;
    double R, I;

    if (type == LAPLACE) f = TWOPI * f;

    for (i = 0; i < 2*nz; i += 2) {
        R = -ze[i];
        I = -ze[i+1];
        if (R == 0.0 && I == 0.0) { mod *= f; pha += PI/2.0; }
        else {
            mod *= sqrt((f+I)*(f+I) + R*R);
            pha += atan2(f+I, R);
        }
    }
    for (i = 0; i < 2*np; i += 2) {
        R = -po[i];
        I = -po[i+1];
        if (R == 0.0 && I == 0.0) { mod /= f; pha -= PI/2.0; }
        else {
            mod /= sqrt((f+I)*(f+I) + R*R);
            pha -= atan2(f+I, R);
        }
    }
    out[0] = mod * cos(pha) * h0;
    out[1] = mod * sin(pha) * h0;
}



/*==================================================================
 *                Response of symetrical FIR filters
 *=================================================================*/
fir_sym_trans(a, na, h0, wsint, out, type)
double *a, h0, wsint, *out;
int na;
int type;
{
    int k;
    double R = 0.0;

    if (type == FIR_SYM_1) {
        for (k=1; k<na; k++) R += a[k] * cos(wsint * k);
        out[0] = (a[0] + 2.0*R) * h0;
        out[1] = 0.;
    }
    else if (type == FIR_SYM_2) {
        for (k=0; k<na; k++) R += a[k] * cos(wsint * (k+0.5));
        out[0] = 2.0*R * h0;
        out[1] = 0.;
    }
}



/*==================================================================
 *                Response of asymetrical FIR filters
 *=================================================================*/
fir_asym_trans(a, na, h0, w, sint, out)
double *a, h0, w, sint, *out;
int na;
{
int k;
double R = 0.0, I = 0.0;
double wsint, y;
double mod, pha;

    wsint = w *sint;
    for (k = 1; k < na; k++) {
        if (a[k] != a[0]) break;
    }
    if (k == na) {
        if (wsint == 0.0) out[0] = 1.;
        else out[0] = (sin(wsint/2.*na) / sin(wsint/2.)) * a[0];
        out[1] = 0;
        return;
    }
    for (k = 0; k < na; k++) {
            y = wsint * k;
            R += a[k] * cos(y);
            I += a[k] * -sin(y);
    }
    mod = sqrt(R*R + I*I);
    pha = atan2(I,R) + (w*(double)((na-1)/2.0)*sint);
    R = mod * cos(pha);
    I = mod * sin(pha);
    out[0] = R * h0;
    out[1] = I * h0;
}


/*==================================================================
 *                Response of IIR filters
 *=================================================================*/
iir_pz_trans(ze, nz, po, np, h0, wsint, out)
double *ze, *po, h0, wsint, *out;
int nz, np;
{
    int i;
    double mod = 1.0, pha = 0.0;
    double r, j1, j2;
    double R, I;
    double c, s;

    c = cos(wsint);
    s = sin(wsint);
    for (i = 0; i < 2*nz; i += 2) {
        R = c + ze[i];
        I = s + ze[i+1];
        mod *= sqrt(R*R + I*I);
        if (R == 0.0 && I == 0.0) pha += 0.0;
        else                      pha += atan2(I, R);
    }
    for (i = 0; i < 2*np; i += 2) {
        R = c + po[i];
        I = s + po[i+1];
        mod /= sqrt(R*R +I*I);
        if (R == 0.0 && I == 0.0) pha += 0.0;
        else                      pha -= atan2(I, R);
    }
    out[0] = mod * cos(pha) * h0;
    out[1] = mod * sin(pha) * h0;
}



/*==================================================================
 *                Complex multiplication
 *=================================================================*/
zmul(a, b, c, d, rp, ip)
double a, b, c, d, *rp, *ip;
{
double r, i;
    r = a*c - b*d;
    i = b*c + a*d;
    *rp = r;
    *ip = i;
}

