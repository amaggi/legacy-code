/* poly.c -- implements polygons and simplices
   see poly.h for further details

   Author  : Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center
*/

#include <stdio.h>
#include <stdlib.h>
#include "set.h"
#include "mem.h"
#include "poly.h"
#include "qhull.h"
#include "geom.h"
#include "globals.h"

/* globals and constants */

int vertex_id= 0;             /* next available id for new vertex */
int ridge_id= 0;              /* next available id for new ridge */
facetT *TAILfacet;         /* sentinel for last facet in facet_list */

/* functions in alphabetical order */

/*-------------------------------------------------
-checkorientation- checks the orientation of all facets
*/
void checkorientation(facetT *facetlist, int fault) {
  facetT *facet;
  realT dist;
  boolT waserror= False;
  
  FORALLfacet_(facetlist) {
    if ((dist= distplane(interior_point, facet)) > -DISTround) {
      if (fault == DATAfault) {
	fprintf(ferr, "qhull error #1: singular input data\n");
	errexit(NULL, NULL, NULL, NULL);
      }
      fprintf(ferr, "qhull error #5: facet f%d has incorrect orientation. Distance to center point %g\n", facet->id, dist);
      errprint(facet, NULL, NULL, NULL);
      waserror= True;
    }
  } 
  if (waserror)
    errexit(NULL, NULL, NULL, NULL);
} /* checkorientation */


/*-------------------------------------------------
-checkplanes- check that each ridge is convex
returns:
  num_nonconvex updated
  num_distplane preserved
*/
void checkplanes(facetT *facetlist, int fault) {
  facetT *facet, *neighbor, **neighborp;
  vertexT *vertex, **vertexp;
  realT dist;
  
  FORALLfacet_(facetlist) {
    FOREACHneighbor_(facet) {
      FOREACHvertex_(neighbor->vertices)
	vertex->seen= False;
    }
    FOREACHvertex_(facet->vertices)
      vertex->seen= True;
    FOREACHneighbor_(facet) {
      FOREACHvertex_(neighbor->vertices) {
	if (!vertex->seen) {
	  vertex->seen= True;
	  if ((dist= distplane(vertex->point, facet)) > -DISTround) {
	    if (fault == DATAfault) {
	      fprintf(ferr, "qhull error #1: singular input data\n");
	      errexit(NULL, NULL, NULL, NULL);
	    }
	    if (IStracing || dist > DISTround) {
	      fprintf(ferr, "qhull error #6: non-convex ridge: v%d of f%d is %g above f%d\n", vertex->id, neighbor->id, dist, facet->id);
	      errprint(facet, NULL, vertex, NULL);
	    }
	    num_nonconvex++;
	  }
	}
      }
    }
  }
} /* checkplanes */


/*-------------------------------------------------
-checkpoints- checks that all points are inside all facets
*/
void checkpoints(facetT *facetlist, pointT *points, int numpoints) {
  boolT errors= False;
  facetT *facet;
  pointT *point, *pointend;
  realT dist;
  
  FORALLpoint_(points, numpoints) {
    FORALLfacet_(facetlist) {
      if ((dist= distplane(point, facet)) > (DELTAvalue + DISTround)) {
	errors= True;
	fprintf(ferr, "qhull error #7: p%d outside f%d, distance = %6.16g\n", pointid_(point), facet->id, dist);
	printpoint(ferr, point);
	printfacet(ferr, facet);
      }
    }
  }
  if (errors) {
    printsummary(ferr, facetlist);
    errexit(NULL, NULL, NULL, NULL);
  }
} /* checkpoints */


/*-------------------------------------------------
-checkpolygon- checks the correctness of the structure
*/
void checkpolygon(facetT *facetlist) {
  facetT *facet, *neighbor, *vertexfacet, **vertexfacetp;
  ridgeT *ridge, **ridgep, *vertexridge, **vertexridgep;
  vertexT *vertex, **vertexp;
  int numvertices, numridges, numfacets;
  boolT waserror= False;
  
  FORALLfacet_(facetlist) {
    numvertices= setsize (facet->vertices);
    if (hull_dim != numvertices) {
      fprintf (ferr, "qhull error #8: %d-d facet f%d has %d vertices\n",
	       hull_dim, facet->id, numvertices);
      errprint(facet, NULL, NULL, NULL);
      waserror= True;
    }
    numfacets= setsize(facet->neighbors);
    if (numfacets != hull_dim) {
      fprintf(ferr, "qhull error #9: %d-d facet f%d has %d neighbors\n",
	      hull_dim, facet->id, numfacets);
      errprint(facet, NULL, NULL, NULL);
      waserror= True;
    }
    numridges= setsize(facet->ridges);
    if (numridges && numridges != hull_dim) {
      fprintf(ferr, "qhull error #10: %d-d facet f%d has %d ridges\n",
	      hull_dim, facet->id, numridges);
      errprint(facet, NULL, NULL, NULL);
      waserror= True;
    }
    FOREACHridge_(facet->ridges) {
      numvertices= setsize (ridge->vertices);
      if (hull_dim - 1 != numvertices) {
	fprintf (ferr, "qhull error #11: %d-d ridge r%d has %d vertices\n",
		 hull_dim, ridge->id, numvertices);
	errprint(NULL, ridge, NULL, NULL);
	waserror= True;
      }
      if (ridge->top != facet && ridge->bottom != facet) { 
	fprintf(ferr, "qhull error #12: facet f%d contains a ridge r%d, but neither top nor bottom facet of r%d is f%d\n", facet->id, ridge->id, ridge->id,
		facet->id);
	errprint(facet, ridge, NULL, NULL);
	waserror= True;
      }
    }
  }
  if (waserror) {
    printsummary(ferr, facetlist);
    errexit(NULL, NULL, NULL, NULL);
  }
} /* checkpolygon */


/*----------------------------------------
-createfacet- creates an oriented facet from a vertex set
returns:
  facet->vertices= vertices
vertices are in inverse processing order
*/
facetT *createfacet(setT *vertices, boolT toporient) {
  facetT *facet;
  int k;
  
  facet= newfacet();
  facet->vertices= vertices;
  facet->toporient= toporient;
  for(k= 0; k < hull_dim; k++)
    updateneighbors(facet, vertices, 0, k);
  trace2((ferr, "createfacet: facet f%d\n", facet->id));
  return (facet);
} /* createfacet */


/*----------------------------------------
-createnewfacet- creates an oriented facet from a vertex set
returns:
  facet->vertices= vertices
*/
facetT *createnewfacet(setT *vertices, boolT toporient) {
  facetT *facet;
  int k;
  
  facet= newfacet();
  facet->vertices= vertices;
  facet->toporient= toporient;
  for(k= 1; k < hull_dim; k++)
    updateneighbors(facet, vertices, 1, k);
  trace2((ferr, "createnewfacet: facet f%d\n", facet->id));
  return (facet);
} /* createnewfacet */


/*----------------------------------------
-createsimplex- creates a simplex from a set of vertices
returns:
   its facet list
*/
facetT *createsimplex(setT *vertices) {
  facetT *facetlist, *facet;
  vertexT *vertex, **vertexp;
  boolT toporient= True;
  
  TAILfacet= newfacet();
  TAILfacet->previous= NULL;
  TAILfacet->next= NULL;
  facetlist= TAILfacet;
  FOREACHvertex_(vertices) {
    facet= createfacet(setnewdel(vertices, vertex), toporient);
    toporient ^= True;
    if (TAILfacet->previous) {
      TAILfacet->previous->next= facet;
      facet->previous= TAILfacet->previous;
      facet->next= TAILfacet;
      TAILfacet->previous= facet;
    }else {
      facetlist= facet;
      facet->next= TAILfacet;
      facet->previous= NULL;
      TAILfacet->previous= facet;
    }
  }
  trace1((ferr, "createsimplex: first facet f%d\n", facet->id));
  if (num_hashitems != 0) {
    fprintf(ferr, "qhull error #13: hash table not empty at the end\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  return (facetlist);
} /* createsimplex */


/*-------------------------------------------------
-deleteinterior- delete interior facets and related structures
returns:
    deletes their outsidesets, interior ridges, interior vertices, normals
*/
void deleteinterior(facetT **facetlist, setT *interior, setT *horizon) {
  facetT *facet, **facetp, *neighbor, **neighborp;
  vertexT *vertex, **vertexp, *nextvertex;
  ridgeT *ridge, **ridgep, *nextridge;
  ridgeT *ridgestofree= NULL;  /* list of ridges to be deleted */
  vertexT *verticestofree= NULL;  /* list of vertices to be deleted */
  
  FOREACHfacet_(interior) {
    FOREACHridge_(facet->ridges)
      ridge->seen= False;
    FOREACHvertex_(facet->vertices) 
      vertex->seen= False;
  }
  FOREACHfacet_(horizon) {
    FOREACHvertex_(facet->vertices)
      vertex->seen= True;
    FOREACHridge_(facet->ridges)
      ridge->seen= True;
  }
  FOREACHfacet_(interior) {
    FOREACHridge_(facet->ridges) {
      if (!ridge->seen) {
	ridge->next= ridgestofree;
	ridgestofree= ridge;
	ridge->seen= True;
      }
    }
    FOREACHvertex_(facet->vertices) { /* vertices to delete are linked using */
      if (!vertex->seen) {            /* vertex->point field */
	vertex->point= (pointT *)verticestofree;
	verticestofree= vertex;
	vertex->seen= True;
      }
    }
    if (facet->previous) {
      facet->previous->next= facet->next;
      facet->next->previous= facet->previous;
    }else {   /* 1st facet in facetlist */
      *facetlist= facet->next;
      facet->next->previous= NULL;
    }
    setfree(&(facet->neighbors));
    setfree(&(facet->ridges));
    setfree(&(facet->vertices));
    setfree(&(facet->outsideset));
    memfree(facet->normal, normal_size);
    memfree(facet, sizeof(facetT));
  }
  for(ridge= ridgestofree; ridge; ridge= nextridge) {
    nextridge= ridge->next;
    setfree(&(ridge->vertices));
    memfree(ridge, sizeof(ridgeT));
  }
  for(vertex= verticestofree; vertex; vertex= nextvertex) {
    nextvertex= (vertexT *)(vertex->point);
    memfree(vertex, sizeof(vertexT));
    num_iteration--;
  }
  trace1((ferr,"deleteinterior:interior facets marked for future deletion\n"));
} /* deleteinterior */


/*-------------------------------------------------
-fliporient- flips the orientation of simplex
assumes that no ridges exist
*/
void fliporient(facetT *facetlist) {
  facetT *facet, *topfacet;
  
  FORALLfacet_(facetlist)
    facet->toporient ^= True;
} /* fliporient */


/*-------------------------------------------------
-makecone- make cone of new facets from point to horizon ridges
returns:
  list of new facets
  all vertex->neighbors= NULL
*/
facetT *makecone(pointT *point, setT *horizon, setT *interior) {
  facetT *facet, **facetp, *facetlist= TAILfacet, *newfacet= NULL;
  facetT *neighbor, **neighborp;
  vertexT *vertex;
  setT *vertices;
  boolT first= True;
  int skip;

  vertex= newvertex(point);
  FOREACHfacet_(horizon) {
    FOREACHneighbor_(facet) {
      if (neighbor->interior == False)
	continue;
      vertices= vertexintersect(facet->vertices, neighbor->vertices, &skip);
      setprepend(&vertices, vertex);
      newfacet=createnewfacet(vertices,(neighbor->toporient)?!(skip%2):skip%2);
      num_newfacets++;
      setreplace(&(facet->neighbors), neighbor, newfacet);
      setappend(&(newfacet->neighbors), facet);
      if (first) {
	facetlist= newfacet;
	first= False;
      }
      TAILfacet->previous->next= newfacet;
      newfacet->previous= TAILfacet->previous;
      newfacet->next= TAILfacet;
      TAILfacet->previous= newfacet;
    }
  } 
  if (num_hashitems != 0) {
    fprintf(ferr, "qhull error #13: hash table not empty at the end\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  trace1((ferr, "makecone: created new facets from point p%d to horizon\n",
	  pointid_(point)));
  return (facetlist);
} /* makecone */


/*----------------------------------------
-newfacet- creates and allocates space for a facet
*/
facetT *newfacet(void) {
  facetT *facet;
  int dummy;
  
  facet= (facetT *)memalloc(sizeof(facetT), &dummy);
  facet->visitid= 0;
  facet->id= facet_id++;
  dummy= hull_dim;
  facet->normal= NULL;
  facet->vertices= setnew(&dummy);
  facet->ridges= NULL;
  facet->neighbors= setnew(&dummy);
  facet->outsideset= NULL;
  facet->furthestdist= 0.00;
  facet->simplicial= True;
  facet->interior= False; 
  facet->horizon= False;
  trace4((ferr, "newfacet: created facet f%d\n", facet->id));
  return (facet);
} /* newfacet */


/*-------------------------------------------------
-newhashtable- returns hash_table,HASHTABLEsize of at least newsize slots
  assumes hash_table is empty on call
*/
void newhashtable(int newsize) {
  int k;
  
  if (hash_table)
    free(hash_table);
  for(HASHTABLEsize= newsize; ; HASHTABLEsize++) {
    for(k= 2; k <= 10; k++)
      if (!(HASHTABLEsize % k))
	break;
    if (k > 10)
      break;
  }
  if(!(hash_table=(hashtableT**)malloc(HASHTABLEsize*sizeof(hashtableT *)))) {
    fprintf(ferr, "qhull error #14: insufficient memory\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  for(k= 0; k < HASHTABLEsize; k++)
    hash_table[k]= NULL;
} /* newhashtable */


/*----------------------------------------
-newridge- creates and allocates space for a ridge
*/
ridgeT *newridge(void) {
  ridgeT *ridge;
  int dummy;

  ridge= (ridgeT *)memalloc(sizeof(ridgeT), &dummy);
  ridge->id= ridge_id++;     
  ridge->seen= True;
  ridge->top= NULL; 
  ridge->bottom= NULL;
  dummy= hull_dim - 1;
  ridge->vertices= setnew(&dummy);
  trace4((ferr, "newridge: created ridge r%d\n", ridge->id));
  return (ridge);
} /* newridge */


/*----------------------------------------
-newvertex- creates and allocates space for a vertex
*/
vertexT *newvertex(pointT *point) {
  vertexT *vertex;
  int dummy;

  num_iteration++;
  vertex= (vertexT *)memalloc(sizeof(vertexT), &dummy);
  vertex->id= vertex_id++;
  vertex->point= point;
  vertex->seen= True;
  trace4((ferr, "newvertex: vertex v%d created\n", vertex->id));
  return (vertex);
} /* newvertex */


/*----------------------------------------
-updateneighbors- links a facet to its neighbors
*/
void updateneighbors(facetT*facet,setT *vertices,int firstindex,int skipindex){
  hashtableT *preventry= NULL, *entry, *newentry;
  int index= 0;
  unsigned h= 0;
  boolT equal, skipone;
  void **elemAp, **elemBp, **skipp;
  
  elemAp= SETelemaddr_(firstindex, void, vertices);
  h= -(unsigned)(SETelem_(skipindex, void, vertices));
  do
    h += (unsigned)(*elemAp++);
  while(*elemAp);
  h %= HASHTABLEsize;
  for(entry= hash_table[h]; entry; entry= entry->next) {
    if (entry->facet == facet) {
      preventry= entry;
      continue;
    }
    equal= True;
    skipone= False;
    index= 0;
    elemAp= SETelemaddr_(0, void, entry->facet->vertices);
    elemBp= SETelemaddr_(0, void, vertices);
    skipp= SETelemaddr_(skipindex, void, vertices);
    do {
      if(elemBp == skipp)
	continue;
      if(*elemAp++ != *elemBp) {
	if (skipone || *elemAp++ != *elemBp) {
	  equal= False;
	  break;
	}
	skipone= True;
      }
    }while(*(++elemBp));
    if (equal)
      break;
    preventry= entry;
  }
  if (entry) {
    num_hashitems--;
    if (preventry)
      preventry->next= entry->next;
    else
      hash_table[h]= entry->next;
    setappend(&(entry->facet->neighbors), facet);
    setappend(&(facet->neighbors), entry->facet);
    memfree(entry, sizeof(hashtableT));
  }else {
    num_hashitems++;
    newentry= (hashtableT *)memalloc(sizeof(hashtableT), &index);
    newentry->facet= facet;
    newentry->next= hash_table[h];
    hash_table[h]= newentry;
  }
} /* updateneighbors */


/*-------------------------------------------------
-vertexintersect- intersects two vertex sets 
  vertices must be ordered
returns:
  set of common points if the sets differ only in one position
*/
setT *vertexintersect(setT *vertexsetA, setT *vertexsetB, int *intersectindex){
  void **elemAp, **elemBp;
  boolT skipone= False;
  int index= 0, size= hull_dim - 1;
  setT *intersection= setnew(&size);
  
  *intersectindex= -1;
  elemAp= SETelemaddr_(0, void, vertexsetA);
  elemBp= SETelemaddr_(0, void, vertexsetB);
  while(*elemAp && *elemBp) {
    while(*elemAp && ((vertexT *)*elemAp)->id >= ((vertexT *)*elemBp)->id) {
      if (*elemAp  == *elemBp) {
	setappend(&intersection, *elemBp);
	index++;
	if (!*(++elemBp))
	  break;
      }
      elemAp++;
    }
    while(*elemBp&&*elemAp&&((vertexT*)*elemBp)->id>=((vertexT*)*elemAp)->id){
      if (*elemBp == *elemAp) {
	setappend(&intersection, *elemAp);
	if (!*(++elemAp))
	  break;
      }else if (!skipone) {
	*intersectindex= index;
	skipone= True;
      }
      elemBp++;
      index++;
    }
  }
  if (*intersectindex == -1)
    *intersectindex= index;
  return(intersection);
} /* vertexintersect */












