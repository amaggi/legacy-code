/*
* $Id: sac_stats.c,v 1.1 2006/02/24 09:07:57 alessia Exp $
*
* Calculates simple statistics on a seismogram without using depmax/depmin.
* Input: number of points, dt, y(t)
* Output: ymin, ymax, ymean, dy/dt, sigma(y)
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "amutils.h"


void sac_stats(int n, float dt, float *y, float *ymin, float *ymax, 
              float *ymean, float *dydt, float *ysigma)
{
   double ysum;
   int i;
   
   // Sanity check
   if (n < 1) {
     fprintf(stderr,"Npts =  %d in sac_stats, exiting subroutine.\n",n);
     return;
   }
   
   // calculate ymin, ymax and ymean with one pass through the seismogram
   *ymin=*ymax=ysum=y[0];
   for (i=1; i<n ; i++) {
     if (y[i] < *ymin) *ymin = y[i];
     if (y[i] > *ymax) *ymax = y[i];
     ysum += y[i];
   }
   *ymean = (float) (ysum / ((float) n));

   // start block to use numerical recipies fit() to calculate gradient
   // and also calculate standard deviation
   {
     float *tmpx, *tmpy, *tmpsig;
     float a,b,chi2,q,siga,sigb;
     int ndata, mwt;
     float y_rtrend;
     
     // allocate memory for temporary arrays
     tmpx=malloc(n*sizeof(float));
     tmpy=malloc(n*sizeof(float));
     tmpsig=malloc(n*sizeof(float));
     
     // prepare temporary arrays
     ysum=0;
     for(i=0;i<n;i++) {
       tmpx[i]=i*dt;
       tmpy[i]=y[i]-*ymean;
     }
     
     ndata=n;
     mwt=0;
     fit_(tmpx, tmpy, &ndata, tmpsig, &mwt, &a, &b, &siga, &sigb, &chi2, &q);
     *dydt=b;
     
     // calculate standard deviation after trend removal
     ysum=0;
     for(i=0;i<n;i++){
       y_rtrend=(tmpy[i] - (a+b*tmpx[i]));
       //fprintf(stderr,"%f\n",y_rtrend);
       ysum += y_rtrend*y_rtrend;
     }
     *ysigma=(float) sqrt(ysum/((float) n));
     
     // free temporary arrays
     free(tmpx);
     free(tmpy);
     free(tmpsig);
     
   }
 

}
