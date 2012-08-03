/*
* $Id: envelope.c,v 1.1 2005/01/14 02:05:38 alessia Exp $
*
* Calculates the envelope of a time series by computing its hilbert transform
*/

#include <stdlib.h>
#include <math.h>
#include "amutils.h"

void envelope(int n, float *x, float *env)
{
   float *hilb;
   int i;

   hilb=malloc(n*sizeof(float));

   hilbert(n, x, hilb);

   for (i=0;i<n;i++)
     env[i]=sqrt(x[i]*x[i]+hilb[i]*hilb[i]);

   free(hilb);
}
