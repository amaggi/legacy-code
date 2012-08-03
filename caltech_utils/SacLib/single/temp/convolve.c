#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dcomplex.h"


/*******************************************************************/

int convolve(const double *data, int ndata,
	      const double *syn,  int nsyn,
	      double *conv)

{
  int nextpow2();
  void cfft();

  int ncorr,n_max,i;
  dcomplex *cdata,*csyn,*ccorr;

  if (ndata <=0 || nsyn <=0) {
    fprintf(stderr,"No data or syn is read in ...\n");
    return -1;
  }
  
  ncorr = ndata+nsyn-1; // length of correlation time series
  n_max = nextpow2(ncorr); // closest 2 power 

  // dynamic allocate array to avoid defining upper limit of the length
  cdata = (dcomplex *) malloc(n_max * sizeof(dcomplex));
  csyn = (dcomplex *) malloc(n_max * sizeof(dcomplex));
  ccorr = (dcomplex *) malloc(n_max * sizeof(dcomplex));

  // set complex data and syn array
  for (i=0; i<ndata; i++) {cdata[i].re = data[i];  cdata[i].im = 0.;}
  for (i=ndata; i<n_max; i++) cdata[i].re = cdata[i].im = 0.;

  for (i=0; i<nsyn; i++) {csyn[i].re = syn[i]; csyn[i].im = 0.;}
  for (i=nsyn; i<n_max; i++) csyn[i].re = csyn[i].im = 0.;

  // fft both complex data and syn array
  cfft(cdata,n_max,1);
  cfft(csyn,n_max,1);

  // in frequency domain, calculate the fft of correlation
  for (i=0;i<n_max;i++) ccorr[i] = dcmult(cdata[i],csyn[i]);

  // fft back to get the correlation time series
  cfft(ccorr,n_max,-1);
  for (i=0;i<ncorr;i++) conv[i] = ccorr[i].re/n_max;
  
  return ncorr;

}
/*******************************************************************/
int convolve_float(const float *data, int ndata,
	      const float *syn,  int nsyn,
	      float *conv)
{
  double *data_double,*syn_double,*conv_double;
  int nconv,i;

  data_double = (double *) malloc(ndata * sizeof(double));
  syn_double = (double *) malloc(nsyn * sizeof(double));
  nconv = ndata + nsyn - 1;
  conv_double = (double *) malloc( nconv * sizeof(double));


  if (data_double == NULL || syn_double == NULL || conv_double == NULL) {
    fprintf(stderr,"Error in allocating memory for convolution \n");
    return -1;}

  for (i=0;i<ndata;i++) {data_double[i] = (double) data[i];}
  for (i=0;i<nsyn;i++) {syn_double[i] = (double) syn[i];}

  if (convolve(data_double,ndata,syn_double,nsyn,conv_double) != nconv) {
    fprintf(stderr,"Error convolve double precision convolution\n");
    return -1;}
  
  for (i=0;i<nconv;i++) {conv[i] = (float) conv_double[i];}
  return nconv;
}


/******************************************************************/

//fortran wraper

int convolve_(const double *data, int *pndata,
	   const double *syn, int *pnsyn,
	   double *corr)
{
  return convolve(data,*pndata,syn,*pnsyn,corr)+1;

}


