#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "dcomplex.h"

// envelope of time series

int envelope(int numsamp,const double *f,double *fenv)

{
  int hilbert();
  int i;
  double *fhilb;
  
  if(numsamp<=0) {
    fprintf(stderr,"%s\n","No data points read ....");
    return (1);
  }

  fhilb = (double *) malloc(numsamp * sizeof(double));

  if (fhilb == NULL) {
    fprintf(stderr,"Incorrect allocation\n");
    return (2);
  }

  hilbert(numsamp,f,fhilb);

  for( i=0; i<numsamp; i++)       
    fenv[i] = sqrt(f[i]* f[i] + fhilb[i]*fhilb[i]);

  free(fhilb);
  return (0);
}

/***********************************************************/

// hilbert transform of time series

int hilbert(int numsamp,const double *f,double *fhilb)

{

  void cfft();
  int find_length();
  void detrend();
  int i, j,npts_max;
  double bpfilt;
  dcomplex imag;
  dcomplex *x;

  if(numsamp<=0) { 
    fprintf(stderr,"%s\n","No data points read ....");
    return (1);
  }

  npts_max = find_length(numsamp);
  if (npts_max == 0) {
    fprintf(stderr,"Too long input time series \n");
    return (2);
  }
  x = (dcomplex *) malloc(npts_max * sizeof(dcomplex));

  if ( x == NULL ) {
    fprintf(stderr,"Incorrect allocation\n");
    return (3);
  }

  imag.re = 0.;  imag.im = 1.0;

  for( i=0; i<numsamp; i++)   
    {
      x[i].re = f[i]; 
      x[i].im = 0.;
    }
  for(i=numsamp; i<npts_max; i++) 
    x[i].re = x[i].im = 0.;

  cfft( x,npts_max,-1 );

  for (i=1; i<=npts_max/2; ++i)
    { 
      x[i] = dcmult(x[i],imag);
    }

  x[0].re = 0.;
  x[0].im = 0.;

  for( i=(npts_max/2)+1; i<npts_max; i++ ) 
    x[i] = dconj(x[npts_max-i]);

  cfft( x,npts_max,1 );

  for( i=0; i<numsamp; i++)       
    fhilb[i] = x[i].re/npts_max;

  detrend(fhilb, numsamp);
  free(x);
  return (0);
}

/*************************************************************/

// remove the mean and trend of time series

void detrend(double *x, int n)

{
  double ds1,ds2,dan,davei,davex,dslope,dai;
  int i, an;
  
  an = n;
  dan=n;
  ds1=0;
  ds2=0;
  
  for(i=0; i<n; ++i)
    {
      ds1 = ds1+ (double) x[i];
      ds2 = ds2 + ds1;
    }
  davei = 0.5 * (dan+1.0);
  davex = ds1/dan;
  dslope = -12.0*(ds2-davei*ds1)/(dan*(dan*dan-1.0));
  for(i=0; i<n; ++i)
    {
      dai = i;
      x[i] = (double) x[i]- davex - dslope*(dai-davei);
    }
  return;
}

/*******************************************************************/

int  find_length(int numsamp)
{
  int npts_max;
  
  if (numsamp > INT_MAX) {
    fprintf(stderr,"Too long time series\n");
    return 0;
  }
  npts_max = 1;
  while( 1 )
    {
      npts_max *= 2;
      if (npts_max >= numsamp) return npts_max;
    }
}

/******************************************************************/

// Below are fortran wrapers


int envelope_(const int *numsamp,const double *f,double *fenv)

{

  return envelope(*numsamp,f,fenv);

}

int hilbert_(const int *numsamp,double *f,double *fhilb)

{

  return hilbert(*numsamp,f,fhilb);

}

