#include <stdio.h>
/*********************************************************************/

/* cross correlation function */

int xcorr(const double *data, int ndata,
	  const double *syn, int nsyn,
	  double *corr, int *nshift)
     // nshift from 1 - ndata to nsyn - 1

{
  int i,j,jmin,jmax,shift,ncorr;
 
  if ( ndata <=0 || nsyn <=0 ) {
    fprintf(stderr,"No data or syn is read in ...\n");
    return 0;
  }
  ncorr = ndata + nsyn -1;
  for ( i=0; i<ncorr; i++ ) { nshift[i] = i - ndata + 1; }

  for ( i=0; i<ncorr; i++) {
    shift = nshift[i];
    jmin = 1 > 1-shift ? 1 : 1-shift; 
    jmax = ndata < nsyn-shift ? ndata : nsyn-shift;
    corr[i] = 0.;
    for (j = jmin - 1; j < jmax; j++) 
      corr[i] += data[j] * syn[j+shift];
  }

  return ncorr;

}
 
// fortran wraper 
int xcorr_(const double *data, int *pndata,
	   const double *syn, int *pnsyn,
	   double *corr, int *nshift)
{
 
  return xcorr(data,*pndata,syn,*pnsyn,corr,nshift);

}
