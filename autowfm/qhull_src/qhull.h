/* qhull.h -- header file for qhull.c that implements Quick Hull
   
   Author  : Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center

   History:
   02/01/93 - started (using code from existing qhull.c) (hh)
   
   Release date:
   05/19/93

   Modified by M. Sambridge for use in fortran interface qhullf (10/4/94)
 
*/


/* ----------------------------------------------
-constants and flags
*/
#define MAXLINE 1000
#define WHITESPACE " \n\t\v\r\f"

/* ======= -functions and procedures- =========== */

/******** functions in alphabetical order ********/

void buildhull(facetT **facetlist, int numpoints);
void errexit(facetT *facet, ridgeT *ridge, vertexT *vertex, pointT *point);
setT *findhorizon(pointT *point, facetT *facet, setT **interior);
facetT *initialhull(setT *vertices);
void initialize(int argc, char *argv[]);
void partitionall(facetT *facetlist, setT *vertices, 
		  pointT *points,int npoints);
void partitionhorizon(facetT *horizonfacet, facetT *newfacets);
void partitionhorizonpoints(facetT *newfacets, setT *horizon);
void partitioninterior(facetT *facetlist, setT *interior);
void partitionpoint(pointT *point, facetT *facetlist);
void qhull(char *inputfile, int *np, int *nd, int *nd_max, int *ndel_max,
           pointT *point, int *ndel, int *v); 
