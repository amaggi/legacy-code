/*
SAC2MT5
	resp2pz
		This routine reads a response file, which may be made up of
		several dated responses (as they come out of rdseed), 
		and collects infromation about the pole-zero response of
		an instrument.  The information is stored in an array of
		Pz structures, one per response.  The routine is passed
		this array and a pointer to the response file.  The routine 
		returns the number of responses.

Alessia Maggi	Completed March 1999
*/

#include "sac2mt5.h"

int resp2pz(Pz poze[], FILE *fpin)
{

int n_resp;
int i,j;
char line[LINE];

const char title[] = "CHANNEL RESPONSE DATA";		/* constant strings */
const char norm[] = "A0 normalization factor:";		/* needed to search */
const char zeros[] = "Number of zeroes:";		/* through the file */
const char poles[] = "Number of poles:";
const char sens[] = "Sensitivity:";


/* find the number of responses in the response file */
n_resp = 0;
while(fgets(line,LINE,fpin)!=NULL)
	if(strstr(line,title)!=NULL)
		n_resp++;

	
/* rewind the response file */

rewind(fpin);


/* get all the relevant information */

for(i=0;i<n_resp;i++){
	while(strstr(fgets(line,LINE,fpin),title)==NULL)
		;
	fgets(line,LINE,fpin);
	sscanf(line,"%*s %*s %s",poze[i].sta);
	fgets(line,LINE,fpin);
	sscanf(line,"%*s %*s %s",poze[i].netwk);
	fgets(line,LINE,fpin);
	sscanf(line,"%*s %*s %s",poze[i].channel);
	fgets(line,LINE,fpin);
	sscanf(line,"%*s %*s %*s %d,%d,%*d:%*d:%*d",&poze[i].year,&poze[i].jday);
	while(strstr(fgets(line,LINE,fpin),norm)==NULL)
		;
	sscanf(line,"%*s %*s %*s %*s %f",&poze[i].A0);  
	fgets(line,LINE,fpin);
	fgets(line,LINE,fpin);
	sscanf(line,"%*s %*s %*s %*s %d",&poze[i].n_zeros);  
	fgets(line,LINE,fpin);
	sscanf(line,"%*s %*s %*s %*s %d",&poze[i].n_poles);  
	fgets(line,LINE,fpin);
	fgets(line,LINE,fpin);
  	for(j=0;j<poze[i].n_zeros;j++)
		sscanf(fgets(line,LINE,fpin),"%*s %*d %f %f %*f %*f",&poze[i].zeros[j].re,&poze[i].zeros[j].im);
	fgets(line,LINE,fpin);
	fgets(line,LINE,fpin);
	for(j=0;j<poze[i].n_poles;j++)
		sscanf(fgets(line,LINE,fpin),"%*s %*d %f %f %*f %*f",&poze[i].poles[j].re,&poze[i].poles[j].im);



	while(strstr(fgets(line,LINE,fpin),sens)==NULL)
                                                      ;
                                                      
	sscanf(line,"%*s %*s %f",&poze[i].sensitivity);  

}



/* return the number of polezero files */

return n_resp;

}
