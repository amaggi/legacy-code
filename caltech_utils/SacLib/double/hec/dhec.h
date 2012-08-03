#ifndef _hec_h
#define _hec_h

#ifndef _dcomplex_h
#define _dcomplex_h
#endif

// spectrum analysis
int  nextpow2(int numsamp);

int denvelope(int numsamp,const double *f,double *fenv);

int dhilbert(int numsamp,const double *f,double *fhilb);

void ddetrend(double *x, int n);

//convolution
int dconvolve(const double *data, int ndata,
	     const double *syn,  int nsyn,
	     double *conv);

//cross-correlation
int dxcorr(const double *data, int ndata,
          const double *syn, int nsyn,
          double *corr, int *nshift);

int dcrosscorr(const double *data, int ndata,
              const double *syn,  int nsyn,
              double *corr);
int dautocorr(const double *data, int ndata,
             double *corr);
int dautocorr2(const double *data, int ndata,
	      double *corr);

//fortran wrapper

int dconvolve_c(const double *data, int *ndata,
              const double *syn,  int *nsyn,
	       double *conv);
int DCONVOLVE_C(const double *data, int *ndata,
              const double *syn,  int *nsyn,
	       double *conv);
int dxcorr_c(const double *data, int *pndata,
           const double *syn, int *pnsyn,
	    double *corr, int *nshift);
int DXCORR_C(const double *data, int *pndata,
           const double *syn, int *pnsyn,
	    double *corr, int *nshift);
int dcrosscorr_c(const double *data, int *pndata,
           const double *syn, int *pnsyn,
		double *corr);
int DCROSSCORR_C(const double *data, int *pndata,
           const double *syn, int *pnsyn,
		double *corr);
int dautocorr_c(const double *data, int *pndata,
	       double *corr);
int DAUTOCORR_C(const double *data, int *pndata,
	       double *corr);
int dautocorr2_c(const double *data, int *pndata,
		double *corr);
int DAUTOCORR2_C(const double *data, int *pndata,
		double *corr);

int denvelope_c(const int *numsamp,const double *f,double *fenv);
int DENVELOPE_C(const int *numsamp,const double *f,double *fenv);
int dhilbert_c(const int *numsamp,double *f,double *fhilb);
int DHILBERT_C(const int *numsamp,double *f,double *fhilb);


#endif
