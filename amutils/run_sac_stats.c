#include <stdio.h>
#include <amutils.h>
#include <sacio.h>

#define MAXPTS 1000000

void ReadSACfile(char *filename,float *x,struct SacHeader *header,int maxpts,int *io_error);


int main(int argc, char* argv[]){
 
  struct SacHeader h;
  float x[MAXPTS];
  int io_err;
  
  float ymin, ymax, ymean, dydt, ysigma;
 
  // Usage
  if (argc != 2 ) {
    fprintf(stderr,"Usage : run_sac_stats sac_file\n");
    return;
  }
  
  // Read the SAC file
  ReadSACfile(argv[1],x,&h,MAXPTS,&io_err);
  
  // Call sac_stats
  sac_stats(h.npts, h.delta, x, &ymin, &ymax, &ymean, &dydt, &ysigma);
  
  fprintf(stdout,"FILE     : %s\n", argv[1]);
  //fprintf(stdout,"YMIN     : %.6g\n", ymin);
  //fprintf(stdout,"YMAX     : %.6g\n", ymax);
  fprintf(stdout,"YMEAN    : %.6g\n", ymean);
  // dy/dt = trend
  fprintf(stdout,"DY/DT    : %.6g\n", dydt);
  // ysigma = sigma(y) after detrend
  fprintf(stdout,"SIGMA(Y) : %.6g\n", ysigma);
 
  return 0;

}
