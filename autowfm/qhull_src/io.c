/* io.c - Input/Output routines of qhull application

   Author:   Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center

   07/09/93 - created (hh)

   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "set.h"
#include "poly.h"
#include "qhull.h"
#include "globals.h"
#include "io.h"
#include "geom.h"


/*-------------------------------------------
-errprint- prints out the information of the erroneous object
    any parameter may be NULL
*/
void errprint(facetT *atfacet, ridgeT *atridge, vertexT *atvertex,
	      pointT *atpoint) {
  if (atfacet) {
    fprintf(ferr, "ERRONEOUS FACET:\n");
    printfacet(ferr, atfacet);
  }
  if (atridge) {
    fprintf(ferr, "ERRONEOUS RIDGE:\n");
    printridge(ferr, atridge);
  }
  if (atvertex) {
    fprintf(ferr, "ERRONEOUS VERTEX:\n");
    printvertex(ferr, atvertex);
  }
  if (atpoint) {
    fprintf(ferr, "ERRONEOUS POINT:\n");
    printpoint(ferr, atpoint);
  }
} /* errprint */


/*-----------------------------------------
-printfacet- prints a facet to stdout
*/
void printfacet(FILE *fp, facetT *facet) {
  facetT *neighbor, **neighborp;
  ridgeT *ridge, **ridgep;
  pointT *point, **pointp, *coords;
  vertexT *vertex, **vertexp;
  int k;

  fprintf(fp, "- f%d\n", facet->id);
  fprintf(fp, "    - orientation: %s\n",((facet->toporient)? "top":"bottom"));
  fprintf(fp, "    - normal: ");
  if (coords= facet->normal) {
    for(k= hull_dim; k; k--)
      fprintf(fp, "%6.2f ", *coords++);
  }
  fprintf(fp, "\n");
  fprintf(fp, "    - offset: %6.2f\n", facet->offset);
  if (facet->outsideset) {
    fprintf(fp, "    - outside set:\n");
    FOREACHpoint_(facet->outsideset)
      printpoint(fp, point);
  }
  fprintf(fp, "    - vertices: ");
  if (facet->toporient) {
    FOREACHvertex_(facet->vertices)
      fprintf(fp, "v%d ", pointid_(vertex->point));
  }else {
    FOREACHvertexreverse12_(facet->vertices)
      fprintf(fp, "v%d ", pointid_(vertex->point));
  }
  fprintf(fp, "\n    - neighboring facets: ");
  FOREACHneighbor_(facet)
    fprintf(fp, "f%d  ", neighbor->id);
  fprintf(fp, "\n");
  if (facet->ridges) {
    fprintf(fp, "    - ridges:\n");
    FOREACHridge_(facet->ridges) {
      printridge(fp, ridge);
      fprintf(fp, "\n");
    } 
  }
} /* printfacet */


/*-----------------------------------------
-printfacets- print list of facets to stdout
*/
void printfacets(FILE *fp, facetT *facetlist) {
  facetT *facet;
  vertexT *vertex, **vertexp;
  
  FORALLfacet_(facetlist) {
    FOREACHvertex_(facet->vertices)
      vertex->seen= False;
  }
  fprintf(fp, "Vertices and their coordinates:\n\n");
  FORALLfacet_(facetlist) {
    FOREACHvertex_(facet->vertices) {
      if (!vertex->seen) {
	vertex->seen= True;
	printvertex(fp, vertex);
      }
    }
  }
  FORALLfacet_(facetlist)
    printfacet(fp, facet);
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(ferr, "\nProblems: \n\n");
    fprintf(ferr, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(ferr, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(ferr, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(ferr, "  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
} /* printfacets */


/*-----------------------------------------
-printincidences- print the facets, each identified by which vertices it
 contains
 
Modified by M. Sambridge for use in fortran interface qhullf (10/4/94)
 
*/
void printincidences(FILE *fp, facetT *facetlist, 
                    int *nd_max, int *nf_max, int *nf, int a[]) {
  vertexT *vertex, **vertexp;
  facetT *facet;
  int k;
  int i= 0;
  int j= 0;
  int numfacets= 0;

  FORALLfacet_(facetlist) {
    if (DELAUNAY && facet->normal[hull_dim-1] > 0.0)
      continue;
    numfacets++;
  }
  /* fprintf(fp, "Number of facets: %d\n", numfacets); */

  *nf = numfacets; 
  if (numfacets > *nf_max) {
     fprintf(ferr, "Maximum number of simplices: %d
         \nNumber calculated          : %d
         \nError detected in printincidences 
         \nMaxmimum number of Delaunay simplices is too small
         \nIncrease parameter nf_max and array sizes \n"
                  , *nf_max, numfacets);
     errexit(NULL, NULL, NULL, NULL);
  }

  FORALLfacet_(facetlist) {
    if (DELAUNAY && facet->normal[hull_dim-1] > 0.0)
      continue;
    if (facet->toporient) {
      FOREACHvertex_(facet->vertices) {
/*	fprintf(fp, "%d ", pointid_(vertex->point)); */
        k = (*nd_max+1)*j + i;
        a[k] = pointid_(vertex->point);
        i++;
      }
    }else {
      FOREACHvertexreverse12_(facet->vertices) {
/*	fprintf(fp, "%d ", pointid_(vertex->point)); */
        k = (*nd_max+1)*j + i;
        a[k] = pointid_(vertex->point);
        i++;
      }
    }
/*    fprintf(fp, " i = %d j=%d k = %d \n",i,j,k); */
/*     fprintf(fp, " \n");  */
    i=0;
    j++;
  }
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(ferr, "\nProblems in convex hull calculation (subroutine qhull): \n\n");
    fprintf(ferr, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(ferr, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(ferr, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(ferr, "  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
} /* printin */


/*----------------------------------------
-printoff2- print a 2-d VECT file containing facetvertices
*/
void printoff2(FILE *fp, facetT *facetlist) {
  vertexT *vertex, **vertexp;
  facetT *facet;
  int k, numvertices= 0, numfacets= 0;
  
  FORALLfacet_(facetlist) {
    numfacets++;
    FOREACHvertex_(facet->vertices)
      numvertices++;
  }
  fprintf(fp, "VECT\n");
  fprintf(fp, "%d %d 1\n", numfacets, numvertices);
  for(k= 0; k < numfacets; k++)
    fprintf(fp, "2 ");
  fprintf(fp, "\n");
  fprintf(fp, "1 ");
  for(k= 1; k < numfacets; k++)
    fprintf(fp, "0 ");
  fprintf(fp, "\n");
  FORALLfacet_(facetlist) {
    if (facet->toporient) {
      FOREACHvertex_(facet->vertices)
	fprintf(fp, "%6.16g %6.16g 0.00 ", vertex->point[0], vertex->point[1]);
    }else {
      FOREACHvertexreverse12_(facet->vertices)
	fprintf(fp, "%6.16g %6.16g 0.00 ", vertex->point[0], vertex->point[1]);
    }
  }
  fprintf(fp, "\n1.00 0.00 0.00 1.00\n");
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(ferr, "\nProblems: \n\n");
    fprintf(ferr, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(ferr, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(ferr, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(ferr, "  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
} /* printoff2 */


/*----------------------------------------
-printoff3- print a 3-d OFF file containing all points and facetvertices
*/
void printoff3(FILE *fp, facetT *facetlist, pointT *points, int numpoints) {
  vertexT *vertex, **vertexp;
  facetT *facet;
  pointT *point, *pointend, *coords;
  int k, numfacets= 0, numridges= 0;

  FORALLfacet_(facetlist) {
    numfacets++;
    SETsize_(facet->neighbors, k);
    numridges += k;
  }
  numridges /= 2;
  fprintf(fp, "OFF\n");
  fprintf(fp, "%d %d %d\n", numpoints, numfacets, numridges);
  FORALLpoint_(points, numpoints) {
    coords= point;
    for(k= hull_dim; k; k--)
      fprintf(fp, "%6.16g ", *coords++);
    fprintf(fp, "\n");
  }
  FORALLfacet_(facetlist) {
    fprintf(fp, "%d ", hull_dim);
    if (facet->toporient) {
      FOREACHvertex_(facet->vertices)
	fprintf(fp, "%d ", pointid_(vertex->point));
    }else {
      FOREACHvertexreverse12_(facet->vertices)
	fprintf(fp, "%d ", pointid_(vertex->point));
    }
    fprintf(fp, "\n");
  }
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(ferr, "\nProblems: \n\n");
    fprintf(ferr, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(ferr, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(ferr, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(ferr, "  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
} /* printoff3 */


/*-------------------------------------------
-printoff4- print a 4-d OFF file of all ridges (general position)
  */
void printoff4 (FILE *fp, facetT *facetlist, pointT *points, int numpoints) {
  facetT *facet, *neighbor, **neighborp;
  vertexT *vertex, **vertexp;
  setT *ridges= NULL, *vertices, **verticesp;
  pointT *point, *pointend, *coords;
  int k, numridges= 0;
  
  FORALLfacet_(facetlist) {
    if (DELAUNAY && facet->normal[hull_dim-1] > 0.0)
      continue;
    FOREACHneighbor_(facet) {
      if (DELAUNAY && neighbor->normal[hull_dim-1] > 0.0)
	continue;
      vertices= vertexintersect(facet->vertices, neighbor->vertices, &k);
      numridges++;
      setappend(&ridges, vertices);
      setdel(&(neighbor->neighbors), facet);
    }
  }
  fprintf(fp, "4OFF\n");
  fprintf(fp, "%d %d 1\n\n", numpoints, numridges);
  FORALLpoint_(points, numpoints) {
    coords= point;
    for (k= hull_dim; k; k--)
      fprintf(fp, "%6.16g ", *coords++);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  if (ridges) {
    for(verticesp= (setT **)&(ridges->e[0]);vertices= *verticesp;verticesp++){
      fprintf(fp, "%d ", hull_dim - 1);
      FOREACHvertex_(vertices)
	fprintf(fp, "%d ", pointid_(vertex->point));
      fprintf(fp, "\n");
    }
  }
  setfree(&ridges);
  setfree(&vertices);
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(ferr, "\nProblems: \n\n");
    fprintf(ferr, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(ferr, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(ferr, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(ferr, "  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
} /* printoff4 */


/*----------------------------------------
-printpoint- prints the information in a point
*/
void printpoint(FILE *fp, pointT *point) {
  int k;
  
  fprintf(fp, "point id: p%d\n", pointid_(point));
  for(k= hull_dim; k; k--)
    fprintf(fp, "%6.2g ", *point++);
  fprintf(fp, "\n");
} /* printpoint */


/*----------------------------------------
-printridge- prints the information in a ridge
*/
void printridge(FILE *fp, ridgeT *ridge) {
  vertexT *vertex, **vertexp;
  
  fprintf(fp, "     - r%d\n", ridge->id);
  fprintf(fp, "           vertices: ");
  FOREACHvertex_(ridge->vertices)
    fprintf(fp, "v%u ", pointid_(vertex->point));
  fprintf(fp, "\n");
  if (ridge->top && ridge->bottom)
    fprintf(fp, "           between f%d and f%d\n",
	    ridge->top->id, ridge->bottom->id);
} /* printridge */


/*----------------------------------------
-printsummary- prints the summary about the computation
*/
void printsummary(FILE *fp, facetT *facetlist) {
  facetT *facet;
  vertexT *vertex, **vertexp;
  int numfacets= 0, numvertices= 0;
  realT worstdistance= 0.00, distance;

  FORALLfacet_(facetlist) {
    FOREACHvertex_(facet->vertices)
      vertex->seen= False;
  }
  FORALLfacet_(facetlist) {
    if (DELAUNAY && facet->normal[hull_dim-1] > 0.0)
      continue;
    numfacets++;
    FOREACHvertex_(facet->vertices) {
      if (!(vertex->seen)) {
        numvertices++; 
        if (CHECKstructure || IStracing) {
          distance= distplane(vertex->point, facet);
	  distance= fabs_(distance);
          maximize_(worstdistance, distance);
	} 
        vertex->seen= True;
      }
    }
  }
  fprintf(fp, "\nResult:\n\n");
  fprintf(fp, "  Number of facets: %d\n", numfacets);
  fprintf(fp, "  Number of vertices: %d\n", numvertices); 
  fprintf(fp, "\nStatistics:\n\n");
  if (CHECKstructure || IStracing)
    fprintf(fp, "  Maximum distance of a vertex from hyperplane: %6.2g\n",
            worstdistance);
  fprintf(fp, "  Number of facets created: %d\n", facet_id);
  fprintf(fp, "  Number of distance tests in visibility determination: %d\n",
	  num_visibility);
  fprintf(fp,"  Number of distance tests in partitioning: %d\n",num_distpart);
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(fp, "\nProblems: \n\n");
    fprintf(fp, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(fp, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(fp, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(fp,"  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
  fprintf(fp, "\n\n");
} /* printsummary */


/*-----------------------------------------
-printtriangles- print hull as list of triangles with their neighbors
*/
void printtriangles(FILE *fp, facetT *facetlist) {
  coordT dist;
  facetT *facet;
  facetT *neighbor, **neighborp;
  setT *vertices;
  int dummy;
  
  if (hull_dim != 3) {
    fprintf(ferr, "qhull: printtriangles() implemented only for 3-d\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  FORALLfacet_(facetlist) {
    fprintf(fp, "facetnumber: %d \n", facet->id);
    FOREACHneighbor_(facet) {
      vertices= vertexintersect(facet->vertices, neighbor->vertices, &dummy);
      dist= pointdist((SETelem_(0, vertexT, vertices))->point,
                      (SETelem_(1, vertexT, vertices))->point);
      fprintf(fp, "%6.3g %d", dist, neighbor->id);
    }
    fprintf(fp, "\n");
  }
  if (num_coplanar + num_nonconvex + num_nearlysingular) {
    fprintf(ferr, "\nProblems: \n\n");
    fprintf(ferr, "  Number of problems: %d\n", 
	    num_coplanar + num_nonconvex + num_nearlysingular);
    fprintf(ferr, "  Number of coplanar points: %d\n", num_coplanar);
    fprintf(ferr, "  Number of non-convex ridges: %d\n", num_nonconvex);
    fprintf(ferr, "  Number of nearly singular facets: %d\n", num_nearlysingular);
  }
} /* printtriangles */


/*-------------------------------------------
-printvertex- prints the information in a vertex
*/
void printvertex(FILE *fp, vertexT *vertex) {
  pointT *point;
  int k;
  
  fprintf(fp, "vertex id: v%d\n", pointid_(vertex->point));
  fprintf(fp, "coordinates: ");
  point= vertex->point;
  for(k= hull_dim; k; k--)
    fprintf(fp, "%6.2g ", *point++); 
  fprintf(fp, "\n");
} /* printvertex */


/*-------------------------------------------
-readpoints- read points from fin into all_points
    fin is lines of coordinates, one per vertex, first line number of points
returns:
    number of points, array of point coordinates, dimension
*/
coordT *readpoints(int *numpoints, int *dimension) {
  coordT *points, *coords, paraboloid;
  char line[MAXLINE+1], token[MAXLINE+1], *tokp, chr, *linep= NULL;
  int diminput, tokcount, pointindex= 0;
  
  while((linep= fgets(line, MAXLINE, fin)) && !strcmp(line, "\n"))
    ;
  if (!linep) {
    fprintf(ferr, "qhull error #0: incorrect input file\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  if (!(diminput= atoi(strtok(line,WHITESPACE)))) {
    fprintf(ferr,"qhull error #16: 1. line specifies the dimension\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  *dimension= DELAUNAY ? diminput+1 : diminput;
  while((linep= fgets(line, MAXLINE, fin)) && !strcmp(line, "\n"))
    ;
  if (!linep) {
    fprintf(ferr, "qhull error #0: incorrect input file\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  if (!(*numpoints= atoi(strtok(line,WHITESPACE)))) {
    fprintf(ferr,"qhull error #17: 2. line specifies the number of points\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  if (!(coords=points=(coordT*)malloc(*numpoints**dimension*sizeof(coordT)))){
    fprintf(ferr, "qhull error #14: insufficient memory\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  tokp= token;
  while (linep= fgets(line, MAXLINE, fin)) {
    if (!strcmp(linep, "\n") || (pointindex++ >= *numpoints))
      continue;
    tokcount= 0;
    paraboloid= 0.0;
    while(chr= *linep++) {
      if (isspace(chr)) {
        if (tokp != token) {
          *tokp= '\0';
          *coords= atof(token);
	  if (DELAUNAY)
	    paraboloid += (*coords) * (*coords);
	  coords++;
          tokp= token;
          tokcount++;
        }
      }else
        *tokp++= chr;
    }
    if (tokp != token) {
      *tokp= '\0';
      *coords= atof(token);
      if (DELAUNAY)
	paraboloid += (*coords) * (*coords);
      coords++;
      tokcount++;
    }
    if (DELAUNAY)
      *coords++= paraboloid;
    if (tokcount != diminput) {
      fprintf(ferr, "qhull error #18: point %d contained %d coordinates\n", 
	      pointindex, tokcount);
      errexit(NULL, NULL, NULL, NULL);
    }
  }
  if (pointindex != *numpoints) {
    fprintf(ferr,"qhull error #19: input contained %d points\n", pointindex);
    errexit(NULL, NULL, NULL, NULL);
  }
  trace1((ferr,"readpoints: read in %d %d-dimensional points\n",
	  pointindex, *dimension));
  return(points);
} /* readpoints */
