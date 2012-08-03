/* mem.h - prototypes for memory management functions

   Author  : Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu

   The Geometry Center

   04/16/93 - initial version (hh)

   WARNING: memalloc returns double-word aligned addresses
   
   */


/* ======= -functions and procedures- =========== */

/******** functions in alphabetical order ********/

void constructsizetable(void);
void *memalloc(int insize, int *outsize);
void memfree(void *pointer, int size);
