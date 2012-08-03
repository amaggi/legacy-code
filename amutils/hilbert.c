/*
* $Id: hilbert.c,v 1.2 2005/01/14 02:05:38 alessia Exp $
*
* Calculates the hilbert transform of data in delta and puts it in kappa. 
*/


#include <stdio.h>
#include "amutils.h"

void hilbert(int n, float *delta, float *kappa)
{
	float d1, d2, d3, d4;
	int i1, i2;

	for (i1 = 0; i1 < n; i1++)
	{
		kappa[i1] = 0.;
		for (i2 = 1; i2 < n; i2++)
		{
			d1 = (i1+i2<n)? delta[i1+i2]: 0.;
			d2 = (i1-i2>=0)? delta[i1-i2]: 0.;
			d3 = (i1+i2+1<n)? delta[i1+i2+1]: 0.;
			d4 = (i1-i2-1>=0)? delta[i1-i2-1]: 0.;

			kappa[i1] -= 0.5 * (d1-d2) / i2 + 0.5 * (d3 - d4) / (i2 + 1);
		}
		kappa[i1] /= PI;
	}
}
