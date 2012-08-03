/*
* $Id: amutils.h,v 1.4 2006/02/24 09:07:57 alessia Exp $
*/

#define PI 3.14159265

extern void fit_(float *x, float *y, int *ndata, float *sig, int *mwt, 
            float *a, float *b, float *siga, float *sigb, float *chi2,
            float *q);


void hilbert(int n, float *delta, float *kappa);
void envelope(int n, float *x, float *env);
void maxima(float *x, int nx, int *maxima , int *nindex, float w_level);
void bracket_maximum(float *x, int nx, int imax , float w_level, int *i1, int *i2) ;
void sac_stats(int n, float dt, float *y, float *ymin, float *ymax, 
              float *ymean, float *dydt, float *ysigma);

