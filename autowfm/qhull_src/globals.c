/* globals.c -- contains all the globals of the qhull application

   Author: Hannu Huhdanpaa
   Email : hannu@geom.umn.edu

   The Geometry Center

   History:
   03/20/93 - created
   
   */

#include <stdio.h>
#include "set.h"
#include "poly.h"
#include "globals.h"


/* ======= -globals- ============================ */

/* user flags */
boolT APPROXhull=0;     /* true 'A' if points within DELTAvalue from facet */
realT DELTAvalue= 0.00; /* value of DELTA from APPROXhull */
boolT BESTfurthest=0;   /*true 'b' if partition points above horizonfacets */
boolT BESTnewfacet= 0;  /*true 'f' if partition points to best, new facet */
boolT CHECKstructure= 0; /* true 'c' if extensive check wanted before output */
boolT DELAUNAY=0;        /* true 'd' if computing DELAUNAY triangulation */
int   IStracing;     /* non-zero if tracing execution, 1 traces the most */
boolT PRINTincidences=0; /* prints the vertex -identifiers of each facet */
boolT PRINToff=0;        /* true 'o' if printing an OFF-file */
boolT PRINTtriangles=0;  /* true 't' if printing triangles */

/* initialized constants */

realT DISTround;     /* maximum round off error, computed in maxmin() */
FILE *fin;           /* pointer to input file */
FILE *fout;          /* pointer to output file */
FILE *ferr;          /* pointer to error file */
pointT *first_point; /* first input point */
int hull_dim;        /* dimension of hull, set by readpoints() */
pointT *interior_point; /* center point of the initial simplex*/
realT MINdenom;      /* min. abs. value for a denominator (prevent overflow) */
realT NEARzero;      /* within roundoff of zero for gausselim */
int normal_size;     /* size in bytes for facet normals */

/* global variables */

int first_newfacet;     /* new facet if facet->id > first_newfacet */
int visit_id;           /* unique id for searching facet neighborhoods */
int facet_id;           /* id of next, new facet from newfacet() */

/* global buffers */

pointT *facet_points;   /* dimXdim coord array of points for sethyperplane */
coordT *gm_matrix;      /* dimXdim matrix for gram_schmidt and gausselim */
coordT **gm_row;        /* array of gm_matrix rows */
int HASHTABLEsize;      /* size of hash_table */
hashtableT **hash_table;/* hash table for connecting facets together */

/* statistics */

int num_coplanar=0;          /* number of coplanar points */
int num_hashitems= 0;        /* number of items in the hash table */ 
int num_nearlysingular= 0;   /* num. of nearly singular facets */
int num_nonconvex= 0;        /* number of non-convex ridges */
int num_distpart= 0;         /* number of distplane() calls in partitioning */
int num_visibility= 0;       /* number of distplane() calls in determining the
                                visibility */
int num_iteration= 0;        /* iteration number */
int num_partition= 0;        /* number of partitioned points per iteration */
int num_newfacets= 0;        /* number of visible facets per iteration */

















