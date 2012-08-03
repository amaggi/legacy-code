/* 

        A fortran interface to C program "qhull" by H. Huhdanpaa.  
        "qhull" is an n-dimensional convex hull program with options for 
	calculating an n-dimensional Delaunay tesselation.
 
	qhullf_(int 	*np, 
        	int     *nd, 
		int     *nd_max, 
		int     *nf_max, 
		int     *mode, 
		double  *points, 
		int 	*nf, 
        	int 	*v) 

        Input:
                np               number of points
                nd 	 	 dimension
                nd_max           maximum number of dimensions 
                nf_max           maximum number of facets. 
				 if mode = 0 facet = Delaunay simplex.
				 if mode = 1 facet = set of nd nodes on hull.
                mode             mode = 0 then Delaunay; mode=1 then Convex hull
                points           position vectors of nodes p(i,j) should contain
				 the i th co-ordinate of the j th point. The 
				 algorithm requires at least nd+1 input points. 
				 Maximum dimensions are points(nd_max,np_max).
 
        Output:
                nfacet         	 number of Delaunay simplices found
                v                list of node indices 
				 Maximum dimensions are v(nd_max+1,np_max).
			 	 If mode = 0: Delaunay calculation.
				    v(i,j) contains the i th vertex of j th 
				    Delaunay simplex (triangle for 2-D),
				    (i=1,nd+1;j=1,nf). Each vertex is a node 
				    numbered from 0 contained in points. For 2-D
				    vertices are stored in counterclockwise 
				    direction.
			 	 If mode = 1: Convex hull is calculated.
				    v(i,j) contains the i th vertex of j th 
				    facet on the convex hull (i=1,nd;j=1,nfacet)
 
        Calling procedure:
	
	If the Delaunay tesselation and the convex hull are required then
	call the routine twice with a different output array, vertices, in
	the calling routine. Note: the sequence must be mode 0 then mode 1.
	(All other sequences are not allowed, e.g. mode 1 then mode 0.)

        Comments on arrays:
 
	Note: The values of v(i,j) are nodes indices ranging from 0 to np-1, 
	but array arguments begin from j=1. 

	The number of non-zero entries in vertices is equal to the number
        of Delaunay simplices (mode=0), or facets on the convex hull (mode=1).
        This value is determined by qhull but it must not exceed nf_max
        otherwise the program is terminated with an error message. In the 
	fortran calling routine vertices must be dimensioned 
	vertices(nd_max+1,nf_max), where nd_max >= nd, nf_max >= nf.
 
        Changes to the original code:
 
        All options to the original C program (see source) produced
        output to standard out and reported errors to standard err. If 
        mode = 0 the fortran routine "qhullf" is equivalent to "qhull d i" 
	which now returns the Delaunay tessellation in the array "vertices",
	if mode = 1, the action is equivalent to "qhull i" and a list of.
	facets on the convex hull is returned in vertices.
 
        Input is not read from stdin. A list of nodes is passed to qhullf
        in array "points".
 
        The file consists of MacIntosh specific code and UNIX specific code.
        If you compile this program for the UNIX environment, keep #define UNIX,
        otherwise, delete this line.
 
						Fortran interface by:
                                                M. Sambridge, RSES ANU, 10/4/94

						Original C version by
   						Author  : Hannu Huhdanpaa
   						Internet: hannu@geom.umn.edu
   						The Geometry Center 02/06/93. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "set.h"
#include "poly.h"
#include "geom.h"
#include "qhull.h"
#include "globals.h"
#include "mem.h"
#include "io.h"
#define UNIX
extern void *buffer_begin; /* pointer to the beginning of current buffer */ 
extern int *indextable;    /* table of freelists, indexed by size */
#ifdef UNIX

char prompt[]= "\n
qhull- convex hull in general dimension\n

    An   - approximate convex hull, delta value n\n\
    c    - check structure\n\
    d    - Delaunay triangulation by lifting to a paraboloid\n\
    f    - in partitioning select the facet the point is furthest above\n\
    b    - partition the horizon points also
    i    - print the facets, each identified by which vertices it contains\n\
    g    - output can be used as input to geomview, available for 2-d - 4-d\n\
    Tn   - tracing on at level n\n\n\
";

/*-------------------------------------------------
-qhullf- fortran interface for qhull (quick hull) program (see above) 
*/

qhullf_(int *np, int *nd, int *nd_max, int *nf_max, 
        int *mode, pointT *input_p, int *nf, int *v) {

  static int arc=3;
  static int arc2=2; 
  static char *arv[3]={"qhull","d","i"};
  static char *arv2[2]={"qhull","i"};
  int i,j,k,kk,kkk,kkkk;
  pointT *points = 0;
  fout= stdout;
  ferr= stderr;
  fin= stdin;
  if(*mode == 0){
    if(hull_dim != 0){
       fprintf(stderr," Delaunay mode called more than once or\n") ;
       fprintf(stderr," after Convex hull mode\n") ;
       fprintf(stderr," Error in qhull: program will crash\n") ;
       errexit(NULL, NULL, NULL, NULL);
    }
/*
					Delaunay initialization
*/
    initialize(arc, arv);
/*
					allocate memory for input points
*/
    k = *np*(*nd+1);

    if (!(points=(pointT*)malloc(*np*(*nd+1)*sizeof(pointT)))){
    fprintf(ferr, "qhull error #14: insufficient memory\n");
    errexit(NULL, NULL, NULL, NULL);
    }

/*					Lift points of dimension nd to
     					points on paraboloid of dimension 
					nd+1 for convex hull algorithm.
*/
    kkk = 0;
    kkkk = *nd;
    for(k=0;k<*np;k++){
      i = k * (*nd_max);
      j = i + *nd - 1;
      points[kkkk] = 0;
      for (kk=i;kk<=j;kk++){
          points[kkk] = input_p[kk]; 
          points[kkkk] = points[kkkk] + (input_p[kk]) * (input_p[kk]);
	  kkk++;
      }
      kkkk = kkkk + *nd + 1;
      kkk++;
    }
    hull_dim = *nd+1;
  }
/*
					Convex hull initialization
*/

  else if(*mode == 1){
/*
					allocate memory for input points
*/
    k = *np*(*nd);

    if (!(points=(pointT*)malloc(*np*(*nd)*sizeof(pointT)))){
    fprintf(ferr, "qhull error #14: insufficient memory\n");
    errexit(NULL, NULL, NULL, NULL);
    }

    initialize(arc2, arv2);
    kkk = 0;
    for(k=0;k<*np;k++){
      i = k * (*nd_max);
      j = i + *nd - 1;
      for (kk=i;kk<=j;kk++){
          points[kkk] = input_p[kk]; 
	  kkk++;
	  /* fprintf(fout," k = %d i = %d j = %d kk = %d kkk = %d\n",k,i,j,kk,kkk); */
      }
    }
    hull_dim = *nd;   
  };

/*
					call qhull to do the work
*/

  qhull(NULL, np, nd, nd_max, nf_max, points, nf, v);
 
/*  free (points); */
  return;
} /* main */
 
#endif
 


/*-------------------------------------------------
-qhull- general dimension convex hull
inputfile != NULL if Mac version

  Modified by M. Sambridge for use in fortran interface qhullf (10/4/94)

*/

void qhull(char *inputfile, int *np, int *nd, int *nd_max,
           int *nf_max, pointT *points, int *nf, int *v) {

  pointT *minx, *maxx;
  int k, numpoints;
  setT *maxpoints, *vertices;
  facetT *facetlist, *facet;
  void *buffer, *nextbuffer;
  int i;
  
  if (inputfile) {  /* Mac */
    fout= fopen("qhull.output", "w");
    ferr= fopen("qhull.trace", "w");
    if (!(fin= fopen(inputfile, "r"))) {
      fprintf(ferr, "qhull error #2: input file could not be opened\n");
      fclose(fout);
      fclose(ferr);
      return;
    }
  }

/*					remove original read in of data
					from standard in and use pointer
					of array passed as an argument
 
     first_point= points= readpoints(&numpoints, &hull_dim); 

*/
     first_point= points;
     numpoints = *np;

  if (hull_dim == 1) {
    fprintf(ferr, "qhull error #3: dimension must be > 1\n");
    return;
  }
  normal_size= hull_dim * sizeof(coordT);
  constructsizetable();
  facet_points= (coordT *)memalloc(hull_dim * hull_dim * sizeof(coordT), &k);
  gm_matrix= (coordT *)memalloc(hull_dim * hull_dim * sizeof(coordT), &k);
  gm_row= (coordT **)memalloc(hull_dim * sizeof(coordT *), &k);
  newhashtable(hull_dim * hull_dim * 10);
  maxpoints= maxmin(points, numpoints, hull_dim, &minx, &maxx);
  k= hull_dim + 1;
  vertices= setnew(&k);
  setprepend(&vertices, newvertex(minx));
  setprepend(&vertices, newvertex(maxx));
  initialvertices(&vertices, maxpoints, points, numpoints, hull_dim+1);
  facetlist= initialhull(vertices);
  partitionall(facetlist, vertices, points, numpoints);
  setfree(&maxpoints);
  setfree(&vertices);
  buildhull(&facetlist, numpoints);
  if (CHECKstructure || IStracing) {
    checkpoints(facetlist, points, numpoints);
    checkpolygon(facetlist);
  }else if (num_coplanar || num_nonconvex || num_nearlysingular) {
    checkorientation(facetlist, ALGORITHMfault);
    checkplanes(facetlist, ALGORITHMfault);
    checkpolygon(facetlist);
  }else
    checkorientation(facetlist, ALGORITHMfault);
  if (PRINTincidences) {
    printincidences(fout, facetlist, nd_max, nf_max, nf, v);
    }
  else if (PRINTtriangles) 
    printtriangles(fout, facetlist);
  else if (PRINToff) {
    switch(hull_dim) {
    case 2:
      printoff2(fout, facetlist);
      break;
    case 3:
      printoff3(fout, facetlist, points, numpoints);
      break;
    case 4:
      printoff4(fout, facetlist, points, numpoints);
      break;
    default:
      printfacets(fout, facetlist);
      break;
    }
  }else 
    printsummary(fout, facetlist);
  FORALLfacet_(facetlist)
   setfree(&(facet->outsideset));
  free(points); 
  free(hash_table);
  for(buffer= buffer_begin; buffer; buffer= nextbuffer) {
    nextbuffer= *((void **)buffer);
    free(buffer);
  }
  free(indextable);
  trace1((ferr, "qhull: algorithm completed\n"));
  if (inputfile) { /* Mac */
    fclose(fin);
    fclose(fout);
    fclose(ferr);
  }
} /* main */


/*-------------------------------------------------
-buildhull- constructs a hull by adding points to a simplex
*/
void buildhull(facetT **facetlist, int numpoints) {
  facetT *facet, *newfacet, *newfacets;
  setT *horizon;         /* set of horizon neighbors */
  setT *interior;        /* set of visible facets for furthest */
  pointT *furthest;
  int dummy;

  trace1((ferr, "buildhull> facetlist %d\n", (*facetlist)->id));
  FORALLfacet_(*facetlist) {
    if (!(furthest= SETlast_(facet->outsideset, dummy)) ||
	(APPROXhull == True && facet->furthestdist < DELTAvalue))
      continue;
    first_newfacet= facet_id; 
    SETdellast_(facet->outsideset, dummy);  
    horizon= findhorizon(furthest, facet, &interior);
    newfacets= makecone(furthest, horizon, interior);
    FORALLnewfacet_(newfacets)
      setfacetplane(newfacet);
    partitioninterior(newfacets, interior);
    if (BESTfurthest)
      partitionhorizonpoints(newfacets, horizon);
    deleteinterior(facetlist, interior, horizon);
    setfree(&horizon);
    setfree(&interior);
  }
  trace1((ferr, "buildhull: completed the hull construction\n"));
} /* buildhull */


/*-------------------------------------------
-errexit- return to system after an error
prints useful information
*/
void errexit(facetT *facet, ridgeT *ridge, vertexT *vertex, pointT *point) {
  errprint(facet, ridge, vertex, point);
  exit(1);
} /*errexit */


/*-------------------------------------------------
-findhorizon- find the horizon and interior for a point that is above the facet
returns:
     non-empty set of horizon facets
     interior: visible facets for the point
     sets facet->interior= True, facet->horizon= False for interior facets
     sets facet->horizon= True, facet->interior= False for horizon facets
*/
setT *findhorizon(pointT *point, facetT *facet, setT **interior) {
  facetT *neighbor, **neighborp;
  setT *horizon= NULL;
  coordT dist;
  int interiorcnt= 1, horizoncnt= 0, n;
  
  trace1((ferr, "findhorizon> for point p%d facet f%d\n",pointid_(point),facet->id));
  *interior= NULL;
  setappend(interior, facet);
  facet->interior= True; facet->horizon= False;
  facet->visitid= ++visit_id;
  SETsize_(facet->neighbors, n);
  for(n= 0; facet= SETelem_(n, facetT, *interior); n++) {
    FOREACHneighbor_(facet) {
      if (neighbor->visitid == visit_id)
	continue;
      neighbor->visitid= visit_id;
      num_visibility++;
      if ((dist= distplane(point, neighbor)) > DISTround) {
	setappend(interior, neighbor);
	interiorcnt++;
	neighbor->interior= True; neighbor->horizon= False;
      }else {
	if (dist > -DISTround) {
	  if (IStracing) 
	    fprintf(ferr, "qhull warning #1: distance from p%d to f%d = %6.16g < DISTround = %6.16g\n", pointid_(point), neighbor->id, dist, DISTround);
	  num_coplanar++;
	}
	setappend(&horizon, neighbor);
	horizoncnt++;
	neighbor->interior= False; neighbor->horizon= True;
      }
    }
  }
  if (!horizon) {
    fprintf(ferr, "qhull error #4: empty horizon\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  SETsize_(horizon, n);
  if ((hull_dim - 1) * n > 2 * HASHTABLEsize) {
    newhashtable(3 * HASHTABLEsize);
  }
  trace1((ferr, "findhorizon: %d horizon facets, %d interior facets\n", 
	  horizoncnt, interiorcnt));
  return (horizon);
} /* findhorizon */


/*-------------------------------------------------
-initialhull- constructs the initial hull as a hull_dim simplex of vertices
  returns:
    facetlist for hull
    interior_point= arithmetic center of vertices
    */
facetT *initialhull(setT *vertices) {
  facetT *facet, *facetlist, *firstfacet;
  
  facetlist= createsimplex(vertices);
  interior_point= getcenter(vertices, hull_dim+1);
  firstfacet= facetlist;
  setfacetplane(firstfacet);
  if (distplane(interior_point, firstfacet) > 0)
    fliporient(facetlist);
  FORALLfacet_(facetlist)
    setfacetplane(facet); 
  checkorientation(facetlist, DATAfault);
  checkpolygon(facetlist);
  checkplanes(facetlist, DATAfault);
  trace1((ferr, "initialhull: simplex constructed\n"));
  return (facetlist);
} /* initialhull */


#ifdef UNIX

/*---------------------------------------------
-initialize- read command line arguments,
    recover qhull_command with defaults
*/
void initialize(int argc, char *argv[]) {
  int i;

  APPROXhull= False;
  BESTfurthest= False;
  BESTnewfacet= False;
  CHECKstructure= False;
  DELAUNAY= False;
  PRINToff= False;
  PRINTincidences= False;
  PRINTtriangles= False;
  IStracing= 0;
  if ((argc == 1) && isatty(0)) {
    fprintf(fout, prompt);
    exit(1);
  }
  for (i=1; i<argc; i++) {
    if (argv[i][0] == '-')
      (argv[i])++;
    switch (argv[i][0]) {
    case 'A':
      if (!isdigit(argv[i][1]))
        fprintf(stderr, "qhull warning #2: no delta value given for arg A\n");
      else {
        DELTAvalue= (realT)atof(&argv[i][1]);
        if (DELTAvalue < 0)
          fprintf(stderr, "qhull warning #3: given delta value negative\n");
        else
          APPROXhull= True;
      }
      break;
    case 'b':
      BESTfurthest= True;
      BESTnewfacet= True;
      break;
    case 'c':
      CHECKstructure= True;
      break;
    case 'd':
      DELAUNAY= True;
      break;
    case 'f':
      BESTnewfacet= True;
      break;
    case 'g':
      PRINToff= True;
      break;
    case 'i':
      PRINTincidences= True;
      break;
    case 't':
      PRINTtriangles= True;
      break;
    case 'T':
      IStracing= atoi(&argv[i][1]);
      break;
    default:
      fprintf (stderr, "qhull warning #4: unknown flag %s\n", argv[i]);
      break;
    }
  }
  if (APPROXhull) {
    if (BESTnewfacet == False) {
      BESTfurthest= True;
      BESTnewfacet= True;
    }
  }
} /* initialize */

#endif

/*-------------------------------------------------
-partitionall- partitions all points into the outsidesets for facetlist
   vertices= set of vertices used by facetlist
*/
void partitionall(facetT *facetlist, setT *vertices,
		  pointT *points, int numpoints){
  setT *vertexpoints= NULL; 
  vertexT *vertex, **vertexp;
  pointT *point, *pointend;
  int upperlimit;
  
  FOREACHvertex_(vertices)
    setadd(&vertexpoints, vertex->point);
  upperlimit= setsize(vertexpoints) - 1;
  FORALLpoint_(points, numpoints)
    if (!setin(vertexpoints, point, upperlimit))
      partitionpoint(point, facetlist);
  setfree(&vertexpoints);
  trace1((ferr, "partitionall: partitioned the points into outside sets\n"));
} /* partitionall */


/*-------------------------------------------------
-partitionhorizon- partitions the outsideset of one horizonfacet
*/
void partitionhorizon(facetT *horizonfacet, facetT *newfacets) {
  int dummy;
  pointT *point, **pointp;
  facetT *facet, *bestfacet;
  coordT dist, bestdist, currentdist;
  boolT furthestdeleted= False;
  pointT *furthest= SETlast_(horizonfacet->outsideset, dummy);
  setT *pointstodelete= NULL;
  
  FOREACHpoint_(horizonfacet->outsideset) {
    bestfacet= horizonfacet;
    num_distpart++;
    bestdist= distplane(point, horizonfacet);
    FORALLfacet_(newfacets) {
      num_distpart++;
      if ((dist= distplane(point, facet)) > bestdist) {
	bestfacet= facet;
	bestdist= dist;
      }
    }
    if (bestfacet != horizonfacet) {
      setappend(&pointstodelete, point);
      if (point == furthest)
	furthestdeleted= True;
      if (bestfacet->furthestdist < bestdist) {
	setappend(&(bestfacet->outsideset), point);
	bestfacet->furthestdist= bestdist;
      }else
	setappend2ndlast(&(bestfacet->outsideset), point);
    }
  }
  FOREACHpoint_(pointstodelete)
    setdel(&(horizonfacet->outsideset), point);
  setfree(&pointstodelete);
  if (furthestdeleted && (furthest=SETlast_(horizonfacet->outsideset, dummy))){
    num_distpart++;
    currentdist= distplane(furthest, horizonfacet);
    FOREACHpoint_(horizonfacet->outsideset) {
      num_distpart++;
      if ((dist= distplane(point, horizonfacet)) > currentdist) {
	furthest= point;
	currentdist= dist;
      }
    }
    setdel(&(horizonfacet->outsideset), furthest);
    setappend(&(horizonfacet->outsideset), furthest);
  }
  trace3((ferr, "partitionhorizon: outsideset f%d repartitioned\n", 
	  horizonfacet->id));
} /* partitionhorizon */


/*-------------------------------------------------
-partitionhorizonpoints- partitions the points of horizon facets if needed
*/
void partitionhorizonpoints(facetT *newfacets, setT *horizon){
  facetT *facet, **facetp;;
  
  FOREACHfacet_(horizon) {
    partitionhorizon(facet, newfacets);
  }
  trace2((ferr, "partitionhorizonpoints: partitioned the horizon points\n"));
} /* partitionhorizonpoints */


/*-------------------------------------------------
-partitioninterior- partitions points to the outside sets of facets
*/
void partitioninterior(facetT *facetlist, setT *interior) {
  facetT *facet, **facetp;
  pointT *point, **pointp;
  
  FOREACHfacet_(interior) {
    if (facet->outsideset) {
      FOREACHpoint_(facet->outsideset) {
	partitionpoint(point, facetlist);
      }
    }
  }
  trace1((ferr,"partitioninterior: partitioned interior points into new facets\n"));
} /* partitioninterior */


/*-------------------------------------------------
-partitionpoint- assigns point to a visible facet in facetlist
    BESTnewfacet= assigns point to facet it is furthest above
    otherwise, picks the first one 

    DELTAvalue == 0.00 in non-approximate case
*/
void partitionpoint(pointT *point, facetT *facetlist) {
  coordT dist, bestdist= DELTAvalue;
  facetT *facet, *bestfacet= NULL;
  
  num_partition++;
  FORALLfacet_(facetlist) {
    num_distpart++;
    if ((dist= distplane(point, facet)) > bestdist) {
      bestfacet= facet;
      bestdist= dist;
      if (!BESTnewfacet && bestdist > DISTround)
	break;
    }
  }
  if (bestfacet) {
    if (bestdist < DISTround) {
      if (IStracing)
	fprintf(ferr, "qhull warning #5: distance from the facet < DISTround, distance: %g, DISTround: %g\n", bestdist, DISTround);
      num_coplanar++;
    }
    if (bestfacet->furthestdist < bestdist) {
      setappend(&(bestfacet->outsideset), point);
      bestfacet->furthestdist= bestdist;
    }else {
      setappend2ndlast(&(bestfacet->outsideset), point);
    }
  }
  trace1((ferr, "partitionpoint: point p%d assigned to a facet f%d\n",
	  pointid_(point), getid_(bestfacet)));
} /* partitionpoint */



