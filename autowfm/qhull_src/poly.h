/* poly.h -- header file for poly.c that constructs a simplex

   Author:   Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center
   
   History:
   02/06/93 - changed the fOREACH-loops (bb & hh)
   01/30/93 - correctly working version ready
   01/27/93 - started (hh)
   */

/*-----------------------------------------------
-constants- these constants are used to determine the appropriate error
message in checkorientation() and checkplanes()
*/

#define ALGORITHMfault 0
#define DATAfault 1


/* ----------------------------------------------
-type names and included files, structures defined below
*/

typedef struct vertexT vertexT;
typedef struct ridgeT ridgeT;
typedef struct facetT facetT;
typedef struct hashtableT hashtableT;
typedef enum {False, True} boolT;
#define coordT double   /* if change to float, make local variables realT */
#define realT double   
#define pointT coordT
#define flagT unsigned int

/* ============ -structures- ====================
   see typedefs above
*/


/* ----------------------------------------------
-facetT- specifies a facet
*/

struct facetT {        /* reals must be first */
  coordT  offset;      /* offset of hyperplane */
  coordT  furthestdist;/* distance to furthest point */
  coordT *normal;      /* normal of hyperplane */
  facetT *previous;    /* pointer to the previous facet in the facetlist */
  facetT *next;        /* pointer to the next facet in the facetlist */
  setT   *outsideset;  /* set of points outside this plane
		          if non-empty, last point is furthest */
  setT   *vertices;    /* vertices belonging to this facet */
  setT   *ridges;      /* ridges of the facet */
  setT   *neighbors;   /* neighbors of the facet */
  int     visitid;     /* visit id, needed by buildhull() */
  int     id:24;       /* unique identifier */
  flagT   toporient:1; /* True if facet has top orientation, False otherwise */
  flagT   simplicial:1;/* True if simplicial facet with neighbors */
  flagT   interior:1;  /* True for interior facets */
  flagT   horizon:1;   /* True for horizon facets */
};


/*----------------------------------------------
-ridgeT- specifies a ridge
*/

struct ridgeT {
  setT   *vertices;    /* vertices belonging to this ridge */
  facetT *top;         /* one facet this ridge is part of */
  facetT *bottom;      /* the other facet this ridge is part of */
  ridgeT *next;        /* can be used to link ridges */
  int     id:24;       /* unique identifier */
  flagT   seen:1;      /* required by fliporient to prevent flipping twice */
};

/* ----------------------------------------------
-vertexT- specifies a vertex
*/

struct vertexT {
  pointT *point;       /* pointer to point structure */
  int     id:24;       /* unique identifier */
  flagT   seen:1;      /* used to eliminate counting a vertex more than once */
};

/* ----------------------------------------------
-hashtableT- hash table entry
*/

struct hashtableT {
  facetT     *facet;     /* pointer to a facet */
  hashtableT *next;      /* pointer to the next hash table entry */
};

/* =========== -macros- ========================= */

#define addto_(p, list)   {(p)->next= (list); (list)= (p);}
#define otherfacet_(r, f) (((r)->top == (f)) ? (r)->bottom : (r)->top)
#define getid_(ridge)     ((ridge) ? (ridge)->id : -1)
#define pointid_(point)   ((point - first_point)/hull_dim)
#define fabs_(a) (((a) < 0) ? -(a):(a))
#define maximize_(maxval, val) {if ((maxval) < (val)) (maxval)= (val);}

/* ----------------------------------------------
-FOREACHxxxx- standard for loops
_FORALLxxxx- standard for loops
*/
#define FOREACHpoint_(points) FOREACHsetelement_(pointT, points, point)
#define FORALLpoint_(points, num) \
  for(point= (points), pointend= (points)+hull_dim*(num); \
      point < pointend; point += hull_dim)
#define FOREACHvertex_(vertices) FOREACHsetelement_(vertexT, vertices,vertex)
#define FOREACHvertexreverse12_(vertices) FOREACHsetelementreverse12_(vertexT, vertices, vertex)
#define FOREACHneighbor_(facet) FOREACHsetelement_(facetT, facet->neighbors, neighbor)
#define FOREACHridgevertex_(vertices) FOREACHsetelement_(vertexT, vertices, ridgevertex)
#define FOREACHridge_(ridges) FOREACHsetelement_(ridgeT, ridges, ridge)
#define FOREACHoldridge_(ridges) FOREACHsetelement_(ridgeT, ridges, oldridge)
#define FOREACHfacetridge_(ridges) FOREACHsetelement_(ridgeT, ridges, facetridge)
#define FOREACHvertexridge_(ridges) FOREACHsetelement_(ridgeT, ridges, vertexridge)
#define FORALLfacet_(facetlist) for(facet=(facetlist);facet->next;facet=facet->next)
#define FORALLnewfacet_(facets) for(newfacet=(facets);newfacet->next;newfacet=newfacet->next)
#define FOREACHfacet_(facets) FOREACHsetelement_(facetT, facets, facet)
#define FOREACHnewfacet_(facetlist) for(newfacet=(facetlist);newfacet;newfacet=newfacet->next)
#define FOREACHinteriorfacet_(facets) FOREACHsetelement_(facetT, facets, interiorfacet)
#define FOREACHvertexfacet_(facets) FOREACHsetelement_(facetT, facets, vertexfacet)

/* ======= -functions and procedures- =========== */

/******** functions in alphabetical order ********/

void checkorientation(facetT *facetlist, int fault);
void checkplanes(facetT *facetlist, int fault);
void checkpoints(facetT *facetlist, pointT *points, int numpoints);
void checkpolygon(facetT *facetlist);
facetT *createfacet(setT *vertices, boolT toporient);
facetT *createnewfacet(setT *vertices, boolT toporient);
facetT *createsimplex(setT *vertices);
void deleteinterior(facetT **facetlist, setT *interior, setT *horizon);
void fliporient(facetT *facetlist);
facetT *makecone(pointT *point, setT *horizon, setT *interior);
void newhashtable(int newsize);
facetT *newfacet(void);
ridgeT *newridge(void);
vertexT *newvertex(pointT *point);
void updateneighbors(facetT *facet, setT *vertices, int firstindex, int skipindex);
setT *vertexintersect(setT *vertexsetA, setT *vertexsetB, int *intersectindex);








