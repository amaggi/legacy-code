#include <stdio.h>

#define NMAX 5000
// IMPORTANT : this number must be larger than npts_max in transform.c

int main() 
{  
  double TimeSeries[NMAX];
  double CrossCorrelation[NMAX];
  double HilbertTransform[NMAX];
  double EnvelopeTransform[NMAX];
  int IndexShift[NMAX];
  char filename[120],*fn;
  int npts,npts_corr,i,nhil,nenv;
  FILE *fd;

  // read in file name
  printf("Input sac file name ...\n");
  //  scanf("%s",filename);
  fn = filename; fn = "sacfile";

  // test sac io functions
  npts = read_sac_double(fn,TimeSeries);
  printf("Total Number of Points in sac file: %d\n",npts);
  
  //test correlation function
  npts_corr=xcorr(TimeSeries,npts,TimeSeries,npts,
			      CrossCorrelation,IndexShift);
  printf("Total Number of Points for auto-correlation: %d\n",npts_corr);
  if (npts_corr == 0) { exit (1); }  
  if (npts_corr > NMAX) { exit(2); }

  // test hilbert transform
  nhil = hilbert(npts,TimeSeries,HilbertTransform);
  printf("Done with hilbert TF calculation : %d\n",nhil);
  if (nhil != 0) {exit(1);}
 
  // test envelope function
  nenv = envelope(npts,TimeSeries,EnvelopeTransform);
  printf("Done with envelope calculation : %d\n",nhil);
  if (nenv != 0) {exit(1);}

  // print out result to files
  printf("Print out Result to files\n");

    fd = fopen("testc.corr","w");
  for (i=0;i<npts_corr;i++)
    fprintf(fd,"%d   %g\n",IndexShift[i],CrossCorrelation[i]);
  fclose(fd);
  
  fd = fopen("testc.env","w");
  for (i=0;i<npts;i++)
    fprintf(fd,"%d  %g %g  %g\n",
	    i,TimeSeries[i],HilbertTransform[i],EnvelopeTransform[i]);
  fclose(fd);

  return 0;

}

  
