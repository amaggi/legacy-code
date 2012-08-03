/* geom.h -- header file for geom.c that implements the geometrical part of
   quick hull application

   Author  : Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center

   History:
   02/09/93 - started (using code from existing qhull.c) (hh)
   */


/* ============ -macros- ======================== */

/*-----------------------------------------------
-det2- 2-d determinate
notes:
    det2 from Carling, "Matrix Inversion", Graphics Gems
*/
#define det2_(a1,a2,b1,b2) ((a1)*(b2) - (a2)*(b1))


/*-----------------------------------------------
-dX, dY, dZ- coordinate differences
*/
#define dX(p1,p2)  (points[(p1) * hull_dim] - points[(p2) * hull_dim])
#define dY(p1,p2)  (points[(p1) * hull_dim + 1] - points[(p2) * hull_dim + 1])
#define dZ(p1,p2)  (points[(p1) * hull_dim + 2] - points[(p2) * hull_dim + 2])
#define dW(p1,p2)  (points[(p1) * hull_dim + 3] - points[(p2) * hull_dim + 3])

/*-----------------------------------------------
-traceN(("format\n", vars));  calls printf if IStracing >= N
*/
#define trace1(args) ((IStracing >= 1) && fprintf args)
#define trace2(args) ((IStracing >= 2) && fprintf args)
#define trace3(args) ((IStracing >= 3) && fprintf args)
#define trace4(args) ((IStracing >= 4) && fprintf args)

/* ======= -functions and procedures- =========== */

/******** functions in alphabetical order ********/

coordT *backsubstitute(int sign, int rows, int columns, boolT *zerodiv);
coordT  distplane(pointT *point, facetT *facet);
realT evaluatediagonal(pointT *point, setT *vertices, int dimension);
void gausselim(int *sign, int rows, int columns, boolT *zerodiv);
pointT *getcenter(setT *vertices, int count);
void gram_schmidt(void);
coordT *gramschmidtnormal(void);
void initgramschmidt(pointT *points);
void initializematrix(pointT *points, int rows, int columns);
void initialvertices(setT **vertices, setT *maxpoints, pointT *points,
		     int numpoints, int pointsneeded);
setT *maxmin(pointT *points, int numpoints, int dimension, pointT **minx,
	     pointT **maxx);
coordT *normalize(coordT *normal);
coordT pointdist(pointT *point1, pointT *point2);
void setfacetplane(facetT *newfacets);
pointT *sethyperplane(int toporient, pointT *points, coordT *offset);



