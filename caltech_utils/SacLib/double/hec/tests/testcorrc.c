#include <stdio.h>
#include "drw_file.h"
#include "dhec.h"

#define NMAX 20000
// IMPORTANT : this number must be larger than npts_max in transform.c

int main() 
{  
  double data[NMAX],syn[NMAX];
  double corr[NMAX],hil[NMAX],env[NMAX],eps=1.0e-4;
  int ind[NMAX];
  double corr2[NMAX],conv[NMAX],time,t0,t0_s,dt,dt_s;
  int npts_corr,zero_index,ndata,nsyn,i;
  FILE *fd;
  

  //read in asc files
  printf("Reading sac files\n");
  dread_sacfile("data.sac",&t0,&dt,&ndata,data);
  dread_sacfile("data.sac",&t0_s,&dt_s,&nsyn,syn);
  if (fabs(dt-dt_s) > eps) {
    printf(" Not same sampling intervale\n");
    return(-1);}

  printf ("Number of data: %d\nNumber of syn: %d\n",ndata,nsyn);

  //test xcorr function
  printf("taking brute force xcorr\n");
  npts_corr=dxcorr(data,ndata,syn,nsyn,corr,ind);
  printf("Total Number of Points for x-correlation: %d\n",npts_corr);
  if (npts_corr == 0) { exit (1); }  
  if (npts_corr > NMAX) { exit(2); }

  //test crosscorr function
  printf("Taking cross-correlation function\n");
  zero_index = dcrosscorr(data,ndata,syn,nsyn,corr2);
  if (zero_index <=0) {
    printf("Error x-correlation\n");
    return(-1);}
  printf("Zero index for x-correlation : %d\n",zero_index);
  
  dwrite_ascfile("corr1.ascii",0,dt,npts_corr,corr);
  dwrite_ascfile("corr2.ascii",0,dt,ndata+nsyn-1,corr2);


  // test hil and env
  printf("Taking hilbert and envelope transform \n");
  if (dhilbert(ndata,data,hil) !=0 | denvelope(ndata,data,env) != 0) {
    printf("Error taking hil and env tf\n");
    return(-1);}
  dwrite_ascfile("hil.ascii",t0,dt,ndata,hil);
  dwrite_ascfile("env.ascii",t0,dt,ndata,env);
 
  //test convolve
  printf("Taking convolution\n");
  if (dconvolve(data,ndata,syn,nsyn,conv) <= 0) {
	  printf("Error taking convolve\n");
	  return(-1);
  }
  dwrite_ascfile("conv.ascii",t0,dt,ndata,conv);
  printf("finish writing \n");
  
  return 0;

  
}

  
