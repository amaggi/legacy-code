/* geom.c -- geometric routines of qhull
   
   Author: Hannu Huhdanpaa 
   Email:  hannu@geom.umn.edu

   The Geometry Center

   02/10/93 - orientation by swapping first 2 vertices if bottom oriented (hh)
   02/09/93 - initial version (hh)
   */
   
#include<stdio.h>
#include<math.h>
#include<float.h>
#include "set.h"
#include "mem.h"
#include "poly.h"
#include "geom.h"
#include "globals.h"
#include "qhull.h"
   
/*-------------------------------------------------
-backsubstitute- performs the back substitution (Golub & van Loan 4.4-9)
     for Ax=b where b= [0,...,0,+-1] PA= LU where last row of A= [0,...,0,1]
     1) Ly=Pb == y=b, 2) Ux=y done below, 3) x is the norm.
*/
coordT *backsubstitute(int sign, int rows, int columns, boolT *zerodiv) {
  int i, j;
  coordT *normal= memalloc(normal_size, &i), *point, *ai, *ak;
  
  point= normal + hull_dim - 1;
  *point--= (sign) ? -1.0 : 1.0;
  for(i= rows - 1; i >= 0; i--) {
    *point= 0.0;
    ai= gm_row[i] + i + 1;
    ak= point+1;
    for(j= i+1; j < columns; j++)
      *point -= *ai++ * *ak++;
    if (fabs_(gm_row[i][i]) < NEARzero) {
      memfree(normal, normal_size);
      *zerodiv= True;
      return NULL;
    }
    *point-- /= (gm_row[i])[i];
  }
  return(normal);
} /* backsubstitute */


/*-------------------------------------------
-distplane- get distance from point to facet
returns:
    positive if point is above facet (i.e., outside)
*/
coordT distplane (pointT *point, facetT *facet) {
  coordT *normalp= facet->normal;
  
  switch(hull_dim){
  case 2:
    return(facet->offset + point[0] * normalp[0] + point[1] * normalp[1]);
    break;
  case 3:
    return(facet->offset+point[0]*normalp[0]+point[1]*normalp[1]+point[2]*normalp[2]);
    break;
  case 4:
    return(facet->offset+point[0]*normalp[0]+point[1]*normalp[1]+point[2]*normalp[2]+point[3]*normalp[3]);
    break;
  case 5:
    return(facet->offset+point[0]*normalp[0]+point[1]*normalp[1]+point[2]*normalp[2]+point[3]*normalp[3]+point[4]*normalp[4]);
    break;
  case 6:
    return(facet->offset+point[0]*normalp[0]+point[1]*normalp[1]+point[2]*normalp[2]+point[3]*normalp[3]+point[4]*normalp[4]+point[5]*normalp[5]);
    break;
  case 7:  
    return(facet->offset+point[0]*normalp[0]+point[1]*normalp[1]+point[2]*normalp[2]+point[3]*normalp[3]+point[4]*normalp[4]+point[5]*normalp[5]+point[6]*normalp[6]);
    break;
  case 8:
    return(facet->offset+point[0]*normalp[0]+point[1]*normalp[1]+point[2]*normalp[2]+point[3]*normalp[3]+point[4]*normalp[4]+point[5]*normalp[5]+point[6]*normalp[6]+point[7]*normalp[7]);
    break;
  default:
    {coordT dist, *coordp; int k;
     dist= facet->offset;
     coordp= point;
     for (k= hull_dim; k ; k--)
       dist += *coordp++ * *normalp++;
     return(dist);
     break;
   }
  }
} /* distplane */


/*-------------------------------------------------
-evaluatediagonal- initializes the point matrix, calls gausselim() to perform
the gaussian elimination on it, and then returns the determinant of the result
matrix
*/
realT evaluatediagonal(pointT *point, setT *vertices, int dimension) {
  pointT *first, *firstpoint, *otherpoint, *a;
  int i, j;
  boolT zerodiv= False;
  realT value;
  
  a= gm_matrix;
  for(i= 0; i < dimension; i++, a+= dimension) 
    gm_row[i]= a;
  a= gm_matrix;
  first= firstpoint= (SETelem_(0, vertexT, vertices))->point;
  for(i= 0; i < dimension; i++)
    *a++= *point++ - *first++;
  for(i= 1; i < dimension; i++) {
    first= firstpoint;
    otherpoint= (SETelem_(i, vertexT, vertices))->point;
    for(j= 0; j < dimension; j++)
      *a++= *otherpoint++ - *first++;
  }
  gausselim(&i, dimension, dimension, &zerodiv);
  if (zerodiv)
    return 0.0;
  value= 1.0;
  for(i= 0; i < dimension; i++)
    value *= (gm_row[i])[i];
  return(fabs_(value));
} /* evaluatediagonal */

  
/*-------------------------------------------------
-gausselim- Gaussian elimination with partial pivoting
  coordT data in gm_matrix (row major), indexed by gm_row
  assumes rows <= columns
returns:
  gm_matrix is upper triangular under gm_row (includes row exchanges)
  sets zerodiv and exits, if a NEARzero pivot occurs
  flips sign for each row exchange and each negative diagonal element
*/
void gausselim(int *sign, int rows, int columns, boolT *zerodiv) {
  coordT *ai, *ak, *rowp, *pivotrow;
  coordT n, pivab, pivot;
  int i, j, k, pivoti;
  
  for(k= 0; k < rows; k++) {
    pivot= fabs_((gm_row[k])[k]);
    pivoti= k;
    for(i= k+1; i < rows; i++) {
      if ((pivab= fabs_((gm_row[i])[k])) > pivot) {
	pivot= pivab;
	pivoti= i;
      }
    }
    if (pivot < NEARzero) {
      *zerodiv= True;
      return;
    }
    if (pivoti != k) {
      rowp= gm_row[pivoti]; 
      gm_row[pivoti]= gm_row[k]; 
      gm_row[k]= rowp; 
      *sign ^= 1;
    }
    pivotrow= gm_row[k] + k;
    pivot= *pivotrow++;  /* signed value, and remainder of row */
    for(i= k+1; i < rows; i++) {
      ai= gm_row[i] + k;
      n= (*ai++)/pivot;
      ak= pivotrow;
      for(j= columns - (k+1); j; j--)
	*ai++ -= n * *ak++;
    }
  }
  for(k= 0; k < rows; k++)
    if ((gm_row[k])[k] < 0)
      *sign ^= 1;
} /* gausselim */


/*----------------------------------------------
-getcenter-  gets arithmetic interior_point of an array of points
*/
pointT *getcenter(setT *vertices, int count) {
  int k;
  pointT *center, *point, *coords;
  vertexT *vertex, **vertexp;
  
  center= (pointT *)memalloc(normal_size, &k);
  for (point= center, k= 0; k < hull_dim; point++, k++) {
    *point= 0.0;
    FOREACHvertex_(vertices)
      *point += vertex->point[k];
  }
  for (point= center, k= hull_dim; k; point++, k--)
    *point /= count;
  trace4((ferr, "getcenter: for %d points\n", count));
  return(center);
} /* getcenter */


/*-------------------------------------------------
-gram_schmidt- implements Gram-Schmidt orthogonalization for general dimension
   uses contents of global matrix gm_matrix as input
   output is also in gm_matrix
notes:
   see Golub & van Loan Algorithm 6.2-2
   performs Gram-Schmidt BY ROWS
*/
void gram_schmidt(void) {
  coordT *a, r, *point, *currentrow;
  int i, j, k;
  
  a= gm_matrix;
  for(k= 0; k < hull_dim; k++, a= currentrow + hull_dim) {
    currentrow= a;
    for(r= 0.0, i= hull_dim; i; i--, a++)
      r += *a * *a;
    if (r < MINdenom) {
      fprintf(ferr, "qhull error #1: singular input data\n");
      errexit(NULL, NULL, NULL, NULL);
    }                               /* should be relative to the matrix */
    if ((r= sqrt(r)) < NEARzero) {  /* condition */
      if (IStracing)
	fprintf(ferr, "qhull warning #6: nearly singular facets\n");
      num_nearlysingular++;  /* turns on check structure */
    }
    a= currentrow;
    for(i= hull_dim; i; i--)
      *a++ /= r;
    for(j= k+1; j < hull_dim; j++, a += hull_dim) {
      r= 0.0;
      point= currentrow;
      for(i= hull_dim; i; i--)
	r += *point++ * *a++;
      for(i= hull_dim; i; i--) 
	*(--a) -= *(--point) * r;
    }
  }
} /* gram_schmidt */


/*-------------------------------------------------
-gramschmidtnormal- returns the normal calculated by Gram-Schmidt
This normal is the last row of gm_matrix.
*/
coordT *gramschmidtnormal(void) {
  int k;
  coordT *a, *normal= (coordT *)memalloc(normal_size, &k);
      
  a= gm_matrix + (hull_dim - 1) * hull_dim;
  for(k= 0; k < hull_dim; k++)
    *(normal + k)= *a++;
  return(normal);
} /* gramschmidtnormal */


/*-------------------------------------------------
-initgramschmidt- initializes the matrix for Gram-Schmidt orthogonalization
  uses lastpoint instead of first point for origin of hyperplane basis
  uses firstpoint - interior_point instead of [0,0,...,1] for initial normal 
  a better algorithm is Householder orthogonalization since error is
  independent of matrix condition.  This implementation only produces an
  error on divide by ~0
  really want to determine rank deficiency by singular value decomposition
  */
void initgramschmidt(pointT *points) {
  pointT *lastpoint, *a, *elem, *point;
  int i, k;

  a= gm_matrix;
  lastpoint= points + (hull_dim - 1) * hull_dim;
  elem= points;
  for(i= 1; i < hull_dim; i++) {
    point= lastpoint;
    for(k= hull_dim; k; k--)
      *a++= *elem++ - *point++;
  }
  elem= points;
  point= interior_point;
  for(k= hull_dim; k; k--)
    *a++= *elem++ - *point++;
} /* initgramschmidt */


/*-------------------------------------------------
-initializematrix- initializes gm_matrix/gm_row for Gaussian elimination
*/
void initializematrix(pointT *points, int rows, int columns) {
  coordT *a, *elem, *point;
  int k, j;
  
  a= gm_matrix;
  for(k= 0; k < columns; k++, a += columns)
    gm_row[k]= a;
  a= gm_matrix;
  elem= points+hull_dim;
  for(k= rows; k; k--) {
    point= points;
    for(j= columns; j; j--)
      *a++= *elem++ - *point++;
  }
} /* initializematrix */

  
/*-------------------------------------------------
-initialvertices- creates a non-singular set of initial points
*/
void initialvertices(setT **vertices, setT *maxpoints, pointT *points,
		     int numpoints, int pointsneeded) {
  pointT *point, **pointp, *pointend, *maxpoint;
  vertexT *vertex, **vertexp;
  int k;
  realT maximum, value;
  boolT inuse;
  
  maxpoint = 0; /* MS */
  
  for(k= 2; k < pointsneeded; k++) {
    maximum= 0.0;
    FOREACHpoint_(maxpoints) {
      inuse= False;
      FOREACHvertex_(*vertices)
	if (vertex->point == point) {
	  inuse= True;
	  break;
	}
      if (inuse)
	continue;
      if ((value= evaluatediagonal(point, *vertices, k)) > maximum) {
	maximum= value;
	maxpoint= point;
      }
    }
    if (maximum) 
      setprepend(vertices, newvertex(maxpoint));
    else {
      maximum= 0.0;
      FORALLpoint_(points, numpoints) {
	inuse= False;
	FOREACHvertex_(*vertices)
	  if (vertex->point == point) {
	    inuse= True;
	    break;
	  }
	if (inuse)
	  continue;
	if ((value= evaluatediagonal(point, *vertices, k)) > maximum) {
	  maximum= value;
	  maxpoint= point;
	}
      }
      if (maximum)
	setprepend(vertices, newvertex(maxpoint));
      else {
	fprintf(ferr, "qhull error #1: singular input data\n");
	errexit(NULL, NULL, NULL, NULL);
      }
    }
  }
} /* initialvertices */


/*-------------------------------------------------
-maxmin- collects the maximum and minimum points of input into a set
*/
setT *maxmin(pointT *points, int numpoints, int dimension, 
	     pointT **minx, pointT **maxx) {
  int k;
  realT maxsum= 0.0, maxcoord, maxmaxcoord= 1.0;
  pointT *minimum, *maximum, *point, *pointend;
  setT *set;
  
  k= dimension + 1;
  set= setnew(&k);
  for(k= 0; k < dimension; k++) {
    minimum= maximum= points;
    FORALLpoint_(points, numpoints) {
      if (maximum[k] < point[k])
	maximum= point;
      else if (minimum[k] > point[k])
	minimum= point;
    }
    maxcoord= (fabs_(maximum[k]) > fabs_(minimum[k])) ? 
      fabs_(maximum[k]) : fabs_(minimum[k]);
    if (k == 0) {
      *minx= minimum;
      *maxx= maximum;
    }
    maximize_(maxmaxcoord, maxcoord);
    maxsum += maxcoord;
    setadd(&set, maximum);
    setadd(&set, minimum);
  }
  trace1((ferr, "maxmin: found the maximum and minimum points\n"));
  /* the following two lines calculate roundoff error and the minimum
     absolute value of denominator according to Lemma 3.2-1 of Golub and
     van Loan "Matrix Computation" */
  DISTround= DBL_EPSILON * dimension * maxsum * 1.01;
  MINdenom= FLT_MIN * maxmaxcoord;
  /* calculation of NEARzero is based on error formula 4.4-13 in the above
     mentioned book, authors say n^3 can be ignored and 10 be used in place of
     rho
     */
  NEARzero= 1600 * maxsum * DISTround;
  return(set);
} /* maxmin */


/*-------------------------------------------------
-normalize- does the normalization
*/
coordT *normalize(coordT *normal) {
  coordT *point, n;
  int k;
  
  for(n= 1.0, point=normal, k= hull_dim - 1; k; point++,k--)
    n += *point * *point;
  for(n= sqrt(n), k= hull_dim; k; k--)
    *point-- /= n;
  return(normal);
} /* normalize */


/*-------------------------------------------
-pointdist- distance between points
*/
coordT pointdist(pointT *point1, pointT *point2) {
  coordT dist, diff;
  int k;
  
  dist= 0.0;
  for (k= hull_dim; k ; k--) {
    diff= *point1++ - *point2++;
    dist += diff * diff;
  }
  return(sqrt(dist));
} /* pointdist */


/*-------------------------------------------------
-setfacetplane- sets the hyperplane for newfacet
*/
void setfacetplane(facetT *newfacet) {
  pointT *points, *coords;
  vertexT *vertex, **vertexp;
  int k;
  
  points= facet_points;
  FOREACHvertex_(newfacet->vertices) {
    coords= vertex->point;
    for(k= hull_dim; k; k--)
      *points++= *coords++;
  } 
  if (newfacet->normal)
    memfree(newfacet->normal, normal_size);
  newfacet->normal= sethyperplane(newfacet->toporient, facet_points, 
				  &(newfacet->offset));
} /* setfacetplane */


/*-------------------------------------------------
-sethyperplane- set normalized hyperplane equation
    assumes oriented simplex
returns:
    offset, normal, 
    if degenerate, normalizes to maximum value
    doesn't work for hull_dim < point_dim
    */
pointT *sethyperplane (int toporient, pointT *points, coordT *offset) {
  coordT norm, *a;
  pointT *normal, *firstpoint, *normalpoint;
  ridgeT *ridgep;
  int k, sign= !toporient;
  boolT zerodiv= False;
  
  if (hull_dim == 2) {
    normal= (coordT *)memalloc(normal_size, &k);
    if ((norm= sqrt(dX(1,0)*dX(1,0) + dY(1,0)*dY(1,0))) < MINdenom) {
      if (IStracing)
	fprintf(ferr, "qhull warning #6: nearly singular facets\n");
      num_nearlysingular++;
    }
    if (!toporient)
      norm= -norm;
    normal[0]= dY(1,0)/norm;
    normal[1]= dX(0,1)/norm;
    *offset= -(points[0]*normal[0]+points[1]*normal[1]);
  }else if (hull_dim == 3) {
    normal= (coordT *)memalloc(normal_size, &k);
    normal[0]= det2_(dY(2,0), dZ(2,0), dY(1,0), dZ(1,0));
    normal[1]= det2_(dX(1,0), dZ(1,0), dX(2,0), dZ(2,0));
    normal[2]= det2_(dX(2,0), dY(2,0), dX(1,0), dY(1,0));
    if ((norm= sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2])) <
	MINdenom) {
      if (IStracing)
	fprintf(ferr, "qhull warning #6: nearly singular facets\n");
      num_nearlysingular++;
    }
    if (!toporient)
      norm= -norm;
    normal[0] /= norm;
    normal[1] /= norm;
    normal[2] /= norm;
    *offset=-(points[0]*normal[0]+points[1]*normal[1]+points[2]*normal[2]);
  }else { 
    initializematrix(points, hull_dim-1, hull_dim);
    gausselim(&sign, hull_dim-1, hull_dim, &zerodiv);
    if (!zerodiv) 
      normal= backsubstitute(sign, hull_dim-1, hull_dim, &zerodiv);
    if (zerodiv) { 
      initgramschmidt(points);
      gram_schmidt();
      normal= gramschmidtnormal();
    }else
      normal= normalize(normal);
    firstpoint= points;
    normalpoint= normal;
    *offset= -(*firstpoint++ * *normalpoint++);
    for(k= hull_dim - 1; k; k--)
      *offset -= *firstpoint++ * *normalpoint++;
  }
  if (IStracing >= 2) {
    fprintf (ferr, "sethyperplane: offset %6.2g normal:", *offset);
    for (k=0; k<hull_dim; k++)
      fprintf(ferr, " %6.2g", normal[k]);
    fprintf(ferr, "\n");
  }
  return (normal);
} /* sethyperplane */
