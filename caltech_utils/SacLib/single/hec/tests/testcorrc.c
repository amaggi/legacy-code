#include <stdio.h>
#include "rw_file.h"
#include "hec.h"

#define NMAX 20000
// IMPORTANT : this number must be larger than npts_max in transform.c

int main() 
{  
  float data[NMAX],syn[NMAX];
  float corr[NMAX],hil[NMAX],env[NMAX],eps=1.0e-4;
  int ind[NMAX];
  float corr2[NMAX],conv[NMAX],time,t0,t0_s,dt,dt_s;
  int npts_corr,zero_index,ndata,nsyn,i;
  FILE *fd;
  

  //read in asc files
  printf("Reading sac files\n");
  read_sacfile("data.sac",&t0,&dt,&ndata,data);
  read_sacfile("syn.sac",&t0_s,&dt_s,&nsyn,syn);
  if (fabs(dt-dt_s) > eps) {
    printf(" Not same sampling intervale\n");
    return(-1);}

  printf ("Number of data: %d\nNumber of syn: %d\n",ndata,nsyn);

  //test xcorr function
  printf("taking brute force xcorr\n");
  npts_corr=xcorr(data,ndata,syn,nsyn,corr,ind);
  printf("Total Number of Points for x-correlation: %d\n",npts_corr);
  if (npts_corr == 0) { exit (1); }  
  if (npts_corr > NMAX) { exit(2); }

  //test crosscorr function
  printf("Taking cross-correlation function\n");
  zero_index = crosscorr(data,ndata,syn,nsyn,corr2);
  if (zero_index <=0) {
    printf("Error x-correlation\n");
    return(-1);}
  printf("Zero index for x-correlation : %d\n",zero_index);
  
  write_sacfile("data.sac","corr1.sac",0,npts_corr,corr);
  write_sacfile("data.sac","corr2.sac",0,ndata+nsyn-1,corr2);


  // test hil and env
  printf("Taking hilbert and envelope transform \n");
  if (hilbert(ndata,data,hil) !=0 | envelope(ndata,data,env) != 0) {
    printf("Error taking hil and env tf\n");
    return(-1);}
  write_sacfile("data.sac","hil.sac",t0,ndata,hil);
  write_sacfile("data.sac","env.sac",t0,ndata,env);
 
  //test convolve
  printf("Taking convolution\n");
  if (convolve(data,ndata,syn,nsyn,conv) <= 0) {
	  printf("Error taking convolve\n");
	  return(-1);
  }
  write_sacfile("data.sac", "conv.sac",t0,ndata+nsyn-1,conv);
  printf("finish writing \n");
  
  return 0;

  
}

  
