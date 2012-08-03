#include "dcomplex.h"
#include "dhec.h"

// fortran wraper 

int dconvolve_c(const double *data, int *ndata,
	      const double *syn,  int *nsyn,
	      double *conv)
{
  return dconvolve(data,*ndata,syn,*nsyn,conv);
}
int DCONVOLVE_C(const double *data, int *ndata,
	      const double *syn,  int *nsyn,
	      double *conv)
{
  return dconvolve(data,*ndata,syn,*nsyn,conv);
}


int dconvolve_c_(const double *data, int *ndata,
	      const double *syn,  int *nsyn,
	      double *conv)
{
  return dconvolve(data,*ndata,syn,*nsyn,conv);
}
int DCONVOLVE_C_(const double *data, int *ndata,
	      const double *syn,  int *nsyn,
	      double *conv)
{
  return dconvolve(data,*ndata,syn,*nsyn,conv);
}


//
int dxcorr_c(const double *data, int *pndata,
 	   const double *syn, int *pnsyn,
	   double *corr, int *nshift)
{
  
  return dxcorr(data,*pndata,syn,*pnsyn,corr,nshift);
  
}
int DXCORR_C(const double *data, int *pndata,
 	   const double *syn, int *pnsyn,
	   double *corr, int *nshift)
{
  
  return dxcorr(data,*pndata,syn,*pnsyn,corr,nshift);
  
}
int dxcorr_c_(const double *data, int *pndata,
 	   const double *syn, int *pnsyn,
	   double *corr, int *nshift)
{
  
  return dxcorr(data,*pndata,syn,*pnsyn,corr,nshift);
  
}
int DXCORR_C_(const double *data, int *pndata,
 	   const double *syn, int *pnsyn,
	   double *corr, int *nshift)
{
  
  return dxcorr(data,*pndata,syn,*pnsyn,corr,nshift);
  
}
//
// return the index of zero-shift in fortran
int dcrosscorr_c(const double *data, int *pndata,
	   const double *syn, int *pnsyn,
	   double *corr)
{
  
  return dcrosscorr(data,*pndata,syn,*pnsyn,corr)+1;

}
int DCROSSCORR_C(const double *data, int *pndata,
	   const double *syn, int *pnsyn,
	   double *corr)
{
  
  return dcrosscorr(data,*pndata,syn,*pnsyn,corr)+1;

}
int dcrosscorr_c_(const double *data, int *pndata,
	   const double *syn, int *pnsyn,
	   double *corr)
{
  
  return dcrosscorr(data,*pndata,syn,*pnsyn,corr)+1;

}
int DCROSSCORR_C_(const double *data, int *pndata,
	   const double *syn, int *pnsyn,
	   double *corr)
{
  
  return dcrosscorr(data,*pndata,syn,*pnsyn,corr)+1;

}

//
int dautocorr_c(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr(data,*pndata,corr)+1;

}
int DAUTOCORR_C(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr(data,*pndata,corr)+1;

}
int dautocorr_c_(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr(data,*pndata,corr)+1;

}
int DAUTOCORR_C_(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr(data,*pndata,corr)+1;

}
//
int dautocorr2_c(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr2(data,*pndata,corr)+1;

}
int DAUTOCORR2_C(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr2(data,*pndata,corr)+1;

}
int dautocorr2_c_(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr2(data,*pndata,corr)+1;

}
int DAUTOCORR2_C_(const double *data, int *pndata,
	   double *corr)
{
 
  return dautocorr2(data,*pndata,corr)+1;

}
//
int denvelope_c(const int *numsamp,const double *f,double *fenv)

{

  return denvelope(*numsamp,f,fenv);

}
int DENVELOPE_C(const int *numsamp,const double *f,double *fenv)

{

  return denvelope(*numsamp,f,fenv);

}
int denvelope_c_(const int *numsamp,const double *f,double *fenv)

{

  return denvelope(*numsamp,f,fenv);

}
int DENVELOPE_C_(const int *numsamp,const double *f,double *fenv)

{

  return denvelope(*numsamp,f,fenv);

}

//
int dhilbert_c(const int *numsamp,double *f,double *fhilb)

{

  return dhilbert(*numsamp,f,fhilb);

}
int DHILBERT_C(const int *numsamp,double *f,double *fhilb)

{

  return dhilbert(*numsamp,f,fhilb);

}
int dhilbert_c_(const int *numsamp,double *f,double *fhilb)

{

  return dhilbert(*numsamp,f,fhilb);

}
int DHILBERT_C_(const int *numsamp,double *f,double *fhilb)

{

  return dhilbert(*numsamp,f,fhilb);

}
