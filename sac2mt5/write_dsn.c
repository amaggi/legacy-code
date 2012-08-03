/*
SAC2MT5
	write_header
		Writes the information from an Event structure into the 
		header of the output .dsn-formatted file.
		FORMAT:
			1-2	i2	Year
			3-4	i2	Month
			5-6	12	Day
			7-8	i2	Hour
			9-10	i2	Minute
			11-13	i3	Seconds*10
			14-18	i5	Latitude*100 	(positive North)
			19-24	i6	Longitude*100 	(positive East)
			25-27	i3	Depth, km
			28-29	i2	Magnitude*10
			30-32	i3	Number of Stations in location
		Note:
			Magnitude and Number of Stations can be set to 0.
	write_dsn
		Writes the information from a Sac structure and a pz 
		structure into the output file for each data file.
		FORMAT:

		* Line 1:	Seismogram Header Line:

			1-4	a4	Station code
			5	i1	Instrument type [
					1 = WWSSN 15-100
					2 = WWSSN 30-100
					3 = WWSSN SP
					4 = GDSN LP
					5 = GDSN SP
					7 = BB P	]
			6-8	a3	Component name (spz, lpz, lpn, lpe,
					bbz; replace center letter with v for 
					velocity-type and with a for 
					acceleration-type instrument response)
			9-18	f10	Magnification of instrument
					[ taken as sensitivity / 10^6]
					negative for reversed polarity
			19-30	a12	Event label
			31-32	i2	Hour of initial time [Initial time = 
					time to which the times in the data 
					list are referenced.  For GDSN this is 
					the time of the first point in the 
					list.]
			33-34	i2	Minute of initial time
			35-39	f5	Second of initial time
			40-44	f5	Sampling interval dt (for GDSN only)
			45-51	f7	Station latitude (+ = N) 
					Optional [if not given, the station
					coordinates should be in file
					M5STATIO.DAT]
			52-58	f7	Station longitude (+ = E) 
					Optional
	
		* Line 2:	Number of data points (and instrument constants
				for GDSN):

			1-5	i5	Number of data points
			6-10	i5	Number of zeros in instrument response
					(for GDSN only)
			11-15	i5	Number of poles in instrument response
					(for GDSN only)
			16-25	d10	Calibration constant [A0] (for GDSN
					only)

		* Line 3:	Instrument Zeros and Poles (for GDSN only):
	
				Free real and imaginary zeros (alternate real
				and imaginary parts)
				Free real and imaginary poles (alternate real
				and imagingary parts)

		* Line 4:	Data points:

			For WWSSN:
				Free time and amplitude pairs, time in 
				seconds * 100, amplitude in microns
			For GDSN:
				Free amplitude in microns

	which_instrument
		Determines the instrument type from the network and component 
		data in the Sac structure
	which_component
		Determines the component name from the component data in the 
		Sac structure

	Note: both which_instrument and which_component are rather fiddly. 
	They may have bugs in them so beware !!  The best thing is to check 
	the output file by hand until you are certain that these routines are 
	working properly.

Alessia Maggi 	Completed March 1999
*/

#include "sac2mt5.h"

int write_header(FILE *fp, Event *quake)
{
int year2;

/* 	find a two digit year from a 4 digit year*/
year2 = quake->year % 100;

/*	write header line 	*/  /* add "()" to lat and lon parts after "(int)" : (int)(quake->lat*100), Onur Tan 21Dec00 */

fprintf(fp,"%2d%2d%2d%2d%2d%03d%5d%6d%3d%2d%3d\n",year2,quake->month,quake->day,quake->hour,quake->min,(int)quake->sec*10,(int)(quake->lat*100),(int)(quake->lon*100),quake->dep,0,0);

}

/* ************************************************************* */

int write_dsn(FILE *fp, Sac *data, Pz *poze)
{
char s[LNAMES];
int inst_type;
int i;
int Nstaname;

/* prepare the information for the new format and write it as you go along */

/* first line */ 

/* abbreviate station name to 4 letters maximum :  1-4 a4 */
strncpy(s,data->staname,5);   /* There is a problem (?) with "4", I changed it to "5"  Onur TAN */ 

/*    fprintf(fp, "%4s",s);     This is original with a bug for 3 characters */


 Nstaname = strlen(s);          /*   Onur TAN */

  if(Nstaname == 3)             /* writing STA NAME in good format */
       fprintf(fp,"%3s ",s);   /*  for 3 or 4 characters          */
  else
       fprintf(fp,"%4s",s); 

/* instrument type : 5 i1 */
/*
	1	=	WWSSN 15-100
	2	=	WWSSN 30-100
	3	=	WWSSN SP
	4	=	GDSN LP
	5	=	GDSN SP
	7	=	BB P
*/
if((inst_type = which_instrument(data))==NOTOK){
	fprintf(stderr,"Cannot allocate instrument type.\n");
	return NOTOK;
}
else
	fprintf(fp,"%1d",inst_type);

/* component name : 6-8 a3 */
if(which_component(s,data)==NOTOK){
	fprintf(stderr,"Cannot allocate component type.\n");
	return NOTOK;
}
else
	fprintf(fp,"%3s",s);

/* magnification of instrument : 9-18 f10   defined as sensitivity / 10^6         */
/* If Sensitivity is NULL  you must use add ADDSENS script  Onur TAN */  


      fprintf(fp,"%10.2e",poze->sensitivity / 1e6);


/* event label: 19--30 a12     leave blank   */

     fprintf(fp,"            ");


/* hour of initial time (see note) : 31-32 i2*/
if(data->nzhour == -12345){
	fprintf(stderr,"Reference hour undefined in the sac file.\n");
	return NOTOK;
}
else
	fprintf(fp,"%2d",data->nzhour);

/* minute of initial time (see note) : 33-34 i2 */
if(data->nzmin == -12345){
	fprintf(stderr,"Reference minute undefined in the sac file.\n");
	return NOTOK;
}
else
	fprintf(fp,"%2d",data->nzmin);

/* second of initial time (see note) : 35-39 f5 */
if(data->nzsec == -12345){
	fprintf(stderr,"Reference second undefined in the sac file.\n");
	return NOTOK;
}
else if(data->nzmsec == -12345){
	fprintf(stderr,"Reference millisecond undefined in the sac file.\n");
	return NOTOK;
}
else
	fprintf(fp,"%5.2f",data->nzsec + (data->nzmsec / 10) * 0.01);

/* sampling interval for dt (for GDSN only) : 40-44 f5 --  Add dt for BB-Z, Onur TAN 21Dec00    */
if(inst_type == 4 || inst_type == 5 || inst_type == 7)
	fprintf(fp,"%5.2f",data->delta);
else
	fprintf(fp,"     ");

/* station latitude */
if(data->stla == -12345){
	fprintf(stderr,"Station latitude undefined in the sac file.\n");
	return NOTOK;
}
else
	fprintf(fp,"%7.2f",data->stla);

/* station longitude */
if(data->stlo == -12345){
	fprintf(stderr,"Station longitude undefined in the sac file.\n");
	return NOTOK;
}
else
	fprintf(fp,"%7.2f\n",data->stlo);

/* second line */

/* number of data points : 1-5 i5 */
fprintf(fp,"%5d",data->npts);

/* number of zeros in instrument response (for GDSN only) : 6-10 i5 */
fprintf(fp,"%5d",poze->n_zeros);           

/* number of poles in instrument response (for GDSN only) : 11-15 i5 */
fprintf(fp,"%5d",poze->n_poles);

/* calibration constant (for GDSN only) : 16-25 d10 
	is A0 from the response file
*/
fprintf(fp,"%10.3e\n",poze->A0);

/* third line */

/* instrument zeros and poles (for GDSN only):
free real and imaginary zeros (alternate real and imaginary parts) 
free real and imaginary poles (alternate real and imaginary parts)
*/
for(i=0;i<poze->n_zeros;i++)
	fprintf(fp,"%10.3e %10.3e\n",poze->zeros[i].re, poze->zeros[i].im);

for(i=0;i<poze->n_poles;i++)
	fprintf(fp,"%10.3e %10.3e\n ",poze->poles[i].re, poze->poles[i].im);

/* fourth line */

/* data points:
for WWSSN:
	free time and amplitude pairs, time in seconds * 100, amplitude in microns
for GDSN
	free amplitude in microns
*/

if(inst_type ==1 || inst_type ==2 || inst_type ==3){
/*	WWSSN	*/
	if(data->leven == 0){ 	/* unevenly spaced data - have x and 
				y separately */
		for(i=0;i<data->npts;i++)
			fprintf(fp, "%10.3e %10.3e\n",data->series[i].x * 100, data->series[i].y);
	}
	else {		/* have to construct x values from beg, npts 
			and delta */
		for(i=0;i<data->npts;i++)
			fprintf(fp, "%10.3e %10.3e\n",(data->beg + i * data->delta) * 100, data->series[i].y);
	}
}
else {
/*	GDSN	*/
	for(i=0;i<data->npts;i++)
		fprintf(fp,"%10.3e\n",data->series[i].y);
}

return OK;
}



/* *************************************************************** */

/* function to determine instrument type from data in sac header */
int which_instrument(Sac *data)
{

if (strcmp(data->network,"wwssn")==0 || strcmp(data->network, "WWSSN")==0){
	if(strncmp(data->component,"lh",2)==0 || strncmp(data->component,"LH",2)==0){

		fprintf(stderr,"Wwssn long period record, which one?\n");
		return NOTOK;
		
		
	}
	else if (strncmp(data->component,"sp",2)==0 || strncmp(data->component,"SH",2)==0)
		return 3;
}

else if (strcmp(data->network,"gdsn")==0 || strcmp(data->network, "GDSN")==0){
	if(strncmp(data->component,"lh",2)==0 || strncmp(data->component,"LH",2)==0)
		return 4;
	else if (strncmp(data->component,"sp",2)==0 || strncmp(data->component,"SH",2)==0)
		return 5;
	else if (strncmp(data->component,"bh",2)==0 || strncmp(data->component, "BH",2)==0)
		return 7;	
	else{
		fprintf(stderr,"Gdsn network with unidentifiable component.\n");
		return NOTOK;
	}
}

else{
	fprintf(stderr,"Unrecognized instrument type %s! \n",data->network);
	return NOTOK;
}

}


/* *************************************************************** */
/* function to determine which component is used */
int which_component(char *s, Sac *data)
{
if (strncmp(data->component,"shz",3)==0 || strncmp(data->component, "SHZ",3)==0){
	strcpy(s,"spz");
	return OK;
}
else if (strncmp(data->component,"she",3)==0 || strncmp(data->component, "SHE",3)==0){
	fprintf(stderr,"Short period horizontal components not supported.\n");
	return NOTOK;
}
else if (strncmp(data->component,"shn",3)==0 || strncmp(data->component, "SHN",3)==0){
	fprintf(stderr,"Short period horizontal components not supported.\n");
	return NOTOK;
}
else if (strncmp(data->component,"lhz",3)==0 || strncmp(data->component, "LHZ",3)==0){
	strcpy(s,"lpz");
	return OK;
}
else if (strncmp(data->component,"lhe",3)==0 || strncmp(data->component, "LHE",3)==0){
	strcpy(s,"lpe");
	return OK;
}
else if (strncmp(data->component,"lhn",3)==0 || strncmp(data->component, "LHN",3)==0){
	strcpy(s,"lpn");
	return OK;
}
else if (strncmp(data->component,"bhz",3)==0 || strncmp(data->component, "BHZ",3)==0){
	strcpy(s,"bbz");
	return OK;
}
else if (strncmp(data->component,"bhe",3)==0 || strncmp(data->component, "BHE",3)==0){
	fprintf(stderr,"Broad band horizontal components not supported.\n");
	return NOTOK;
}
else if (strncmp(data->component,"bhn",3)==0 || strncmp(data->component, "BHN",3)==0){
	fprintf(stderr,"Broad band horizontal components not supported.\n");
	return NOTOK;
}
else{
	fprintf(stderr,"Unrecognized component %s.\n",data->component);
	return NOTOK;
}

}
