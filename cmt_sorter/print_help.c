#include"header.h"

void print_help(){

  fprintf(stderr,"\n");
  fprintf(stderr,"Program CMT_sorter - (c) Steve Mangino, Keith Priestley, Alessia Maggi 1994-2002\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Purpose: compile a list of events conforming to search criteria.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Usage: CMT_sorter [parameters]\n");
  fprintf(stderr,"\n");
  fprintf(stderr," -help : Print this help\n");
  fprintf(stderr,"\n");
  fprintf(stderr," -cat [cmt|pde|ehb|txt|jdy] <filename> : file to search with type\n");
  fprintf(stderr," -cen : use CMT hypocenter coordinates\n");
  fprintf(stderr,"\n");
  fprintf(stderr," -mom <minmo> <maxmo> : exponent of seismic moment search\n");
  fprintf(stderr," -mag <minmag> <maxmag> : magnitude search\n");
  fprintf(stderr," -ms  <minmag> <maxmag> : surface wave magnitude search\n");
  fprintf(stderr," -mb  <minmag> <maxmag> : body wave magnitude search\n");
  fprintf(stderr,"\n");
  fprintf(stderr," -dep <mindep> <maxdep> : depth search\n");
  fprintf(stderr," -cir <cenlat> <cenlon> <rmin> <rmax> : circle search\n");
  fprintf(stderr," -seg <cenlat> <cenlon> <rmin> <rmax> <bazmin> <bazmax> : circle segment search\n");
  fprintf(stderr," -rec <minlat> <maxlat> <minlon> <maxlon> : rectangle search\n");
  fprintf(stderr," -date <yyyymmdd> <hhmm> <yyyymmdd> <hhmm> : date search\n");
  fprintf(stderr,"\n");
  fprintf(stderr," -noh : do not include header in .txt file\n");
  fprintf(stderr," -lhd : long header, include command line in .txt file\n");
  fprintf(stderr," -reg : include region in .txt file\n");
  fprintf(stderr," -jdy : use different output to .txt file, including julian day\n");
  fprintf(stderr," -out <filename> : base name of output file, default is search\n");
  fprintf(stderr,"\n");
  fprintf(stderr," N.B.1 : If more than one conflicting flag is specified, only the last one will be considered\n");
  fprintf(stderr,"\n");
  fprintf(stderr," N.B.2 : Paths to catalogs in madingley.org: /data/teleseis/Catalogs/... \n");

}
