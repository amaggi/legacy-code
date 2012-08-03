/*
SAC2MT5:
	read_sac:
		Reads the header and data parts of an alpha formatted sac
		file.  It stores the information in a structure of type
		Sac defined in sac2mt5.h.  The routine is passed a pointer 
		to a Sac structure and a pointer to the sac file.  The
		routine returns NOTOK if there has been a problem reading
		the sac file.
	read_txt:
		Reads event information from a .txt file in the -jday format
		with one BLANK header line.  The information is stored in an 
		Event structure, defined in sac2mt5.h.  The routine is passed 
		a pointer to an Event structure and a pointer to the .txt 
		file.  The routine returns NOTOK if there has been a problem 
		reading the .txt file.
	read_list:
		Reads a file containing a paired list of data and response 
		filenames, stores these filenames in an array of 
		Npair structures (essentially two filenames) and returns the 
		number of such pairs.  The routine is passed this array and a 
		pointer to the list-file, which is formatted as follows:
			data_name1	resp_name1
			data_name2	resp_name2
			..........	.........
		with no blank line at the end of the file.

Alessia Maggi	Completed March 1999.

*/

#include"sac2mt5.h"


int read_sac(Sac *data, FILE *fp)
{
int n_lines;	/* number of lines of data */
int n_last;	/* number of columns of data in the last line */
int i;

/* read header information of interest into the Sac structure data */



fscanf(fp,"%f %*f %*f %*f %*f",&data->delta);
fscanf(fp,"%f %f %f %*f %*f",&data->beg,&data->ed,&data->origin);
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %f %f %*f %*f",&data->stla,&data->stlo);
fscanf(fp,"%f %f %*f %*f %*f",&data->evla,&data->evlo);
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%*f %*f %*f %*f %*f");
fscanf(fp,"%d %d %d %d %d",&data->nzyear,&data->nzjday, &data->nzhour, &data->nzmin, &data->nzsec);
fscanf(fp,"%d %*d %*d %*d %d",&data->nzmsec,&data->npts);
fscanf(fp,"%*d %*d %*d %*d %*d");
fscanf(fp,"%*d %*d %*d %*d %*d");
fscanf(fp,"%*d %*d %*d %*d %*d");
fscanf(fp,"%*d %*d %*d %*d %*d");
fscanf(fp,"%*d %*d %*d %*d %*d");
fscanf(fp,"%d %*d %*d %*d %*d", &data->leven);
fscanf(fp,"%s %*s", &data->staname);
fscanf(fp,"%*s %*s %*s");
fscanf(fp,"%*s %*s %*s");
fscanf(fp,"%*s %*s %*s");
fscanf(fp,"%*s %*s %*s");
fscanf(fp,"%*s %*s %*s");
fscanf(fp,"%*s %*s %s", &data->component);
fscanf(fp,"%s %*s %*s",&data->network);

n_lines = data->npts / 5; 	
n_last = data->npts % 5;	


/*	read the dependent variable 	*/

for(i=0;i<n_lines;i++)		/* read all the complete lines */
{
	fscanf(fp,"%f %f %f %f %f", &data->series[5*i].y, &data->series[5*i+1].y, &data->series[5*i+2].y, &data->series[5*i+3].y, &data->series[5*i+4].y);
}
/********Onur TAN, 21Dec2000


If n_last=0, the swich statement gives error because f the 'default' value.
I add 'case 0' to fix it. 

********/



switch(n_last){		/* read the last line */
        case 0 : 
                break;
        case 1 : 
		fscanf(fp,"%f",&data->series[5*n_lines].y);
		break;
	case 2 : 
		fscanf(fp,"%f %f",&data->series[5*n_lines].y, &data->series[5*n_lines+1].y);
		break;
	case 3 : 
		fscanf(fp,"%f %f %f",&data->series[5*n_lines].y, &data->series[5*n_lines+1].y, &data->series[5*n_lines+2].y);
		break;
	case 4 : 
		fscanf(fp,"%f %f %f %f",&data->series[5*n_lines].y, &data->series[5*n_lines+1].y, &data->series[5*n_lines+2].y, &data->series[5*n_lines+3].y);
		break;
	default : 
		fprintf(stderr,"Error in reading y data in sac file.\n");
		return NOTOK;
		break;
}




if(data->leven == 0 ) {		/* data is not evenly spaced, so both 
				dependent and independent variables are 
				explicitly listed in the sac file  - the 
				independent variable is listed after the 
				dependent one */

for(i=0;i<n_lines;i++)		/* read all the complete lines */
	fscanf(fp,"%f %f %f %f %f", &data->series[5*i].x, &data->series[5*i+1].x, &data->series[5*i+2].x, &data->series[5*i+3].x, &data->series[5*i+4].x);

switch(n_last){			/* read the last line */
	case 1 : fscanf(fp,"%f",&data->series[5*n_lines].x);
	case 2 : fscanf(fp,"%f %f",&data->series[5*n_lines].x, &data->series[5*n_lines+1].x);
	case 3 : fscanf(fp,"%f %f %f",&data->series[5*n_lines].x, &data->series[5*n_lines+1].x, &data->series[5*n_lines+2].x);
	case 4 : fscanf(fp,"%f %f %f %f",&data->series[5*n_lines].x, &data->series[5*n_lines+1].x, &data->series[5*n_lines+2].x, &data->series[5*n_lines+3].x);
	default : {
		fprintf(stderr,"Error in reading x data in sac file.\n");
		return NOTOK;
		break;
	}
}
}

return OK;
}




/* ************************************************************* */
/* read information from a .txt file */

int read_txt(FILE *fp, Event *quake)
{
fscanf(fp,"%d (%d) %d %d %d %d %f %f %f %d %f %*d",&quake->year, &quake->jday, &quake->month, &quake->day, &quake->hour, &quake->min, &quake->sec, &quake->lat, &quake->lon, &quake->dep, &quake->mag);
/*fprintf(stderr,"%d (%d) %d %d %d %d %f %f %f %d %f %*d",&quake->year, &quake->jday, &quake->month, &quake->day, &quake->hour, &quake->min, &quake->sec, &quake->lat, &quake->lon, &quake->dep, &quake->mag); */

return OK;
}


/* ************************************************************* */
/* find the names of data-response pairs */

int find_pairs(FILE *fp, Npair names[])
{
int i, n_names;

i=0;
while(fscanf(fp,"%s %s",names[i].data, names[i].resp) != EOF)
	i++;

n_names = i;
return n_names;

}
