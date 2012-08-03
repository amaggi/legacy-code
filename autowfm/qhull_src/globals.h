/* globals.h -- contains all the globals of the qhull application as externs
   used by other modules of this appication, see globals.c for details
   
   Author: Hannu Huhdanpaa
   Email : hannu@geom.umn.edu

   The Geometry Center

   History:
   03/20/93 - created
   
   */

/* ======= -globals- ============================ */

extern boolT        APPROXhull; 
extern boolT        BESTfurthest;
extern boolT        BESTnewfacet;
extern boolT        CHECKstructure;
extern boolT        DELAUNAY;
extern realT        DELTAvalue;
extern realT        DISTround;
extern int          facet_id;  
extern FILE        *fin;
extern FILE        *fout;
extern FILE        *ferr;
extern pointT      *facet_points;
extern int          first_newfacet;
extern pointT      *first_point;
extern coordT      *gm_matrix; 
extern coordT     **gm_row;
extern int          HASHTABLEsize;
extern hashtableT **hash_table;
extern int          hull_dim;
extern pointT      *interior_point;
extern int          IStracing;
extern realT        MINdenom;
extern realT        NEARzero;
extern int          normal_size;
extern int          num_coplanar;
extern int          num_hashitems;
extern int          num_nearlysingular;
extern int          num_nonconvex;
extern int          num_distpart;
extern int          num_iteration;
extern int          num_partition;
extern int          num_newfacets;
extern int          num_visibility;
extern boolT        PRINTincidences;
extern boolT        PRINToff; 
extern boolT        PRINTtriangles;
extern int          visit_id;

