#include <stdio.h>

#define NMAX 5000
// IMPORTANT : this number must be larger than npts_max in transform.c

int main() 
{  
  double TimeSeries1[NMAX],TimeSeries2[NMAX];
  double CrossCorrelation[NMAX];
  int IndexShift[NMAX];
  double CrossCorrelation2[NMAX];
  int npts_corr,zero_index,ndata,nsyn,i;
  FILE *fd;
  

  fd = fopen("data.tmp","r");
  i = 0;
  while (fscanf(fd,"%lf",&(TimeSeries1[i])) != EOF) i++;
  ndata = i;
  close(fd);
  fd = fopen("syn.tmp","r");
  i = 0;
  while (fscanf(fd,"%lf",&(TimeSeries2[i])) != EOF) i++;
  nsyn = i;
  close(fd);

  printf ("Number of data: %d\nNumber of syn: %d\n",ndata,nsyn);

  //test xcorr function
  npts_corr=xcorr(TimeSeries1,ndata,TimeSeries2,nsyn,
		  CrossCorrelation,IndexShift);
  printf("Total Number of Points for auto-correlation: %d\n",npts_corr);
  if (npts_corr == 0) { exit (1); }  
  if (npts_corr > NMAX) { exit(2); }

  zero_index = crosscorr(TimeSeries1,ndata,TimeSeries2,nsyn,
			 CrossCorrelation2);

  // print out result to files
  printf("Print out Result to files\n");
  fd = fopen("testc.corr1","w");
  for (i=0;i<npts_corr;i++)
    fprintf(fd,"%d   %g\n",IndexShift[i],CrossCorrelation[i]);
  fclose(fd);
  
  fd = fopen("testc.corr2","w");
  for (i=0;i<ndata+nsyn-1;i++) {
    fprintf(fd,"%d    %g\n",i-zero_index,CrossCorrelation2[i]);
  }
  fclose(fd);
  return 0;

}

  
