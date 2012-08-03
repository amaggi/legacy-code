/* mem.c - memory management routines for qhull
   
   Author: Hannu Huhdanpaa
   Internet: hannu@geom.umn.edu
   
   The Geometry Center
   
   04/16/93 - initial version (hh)
   
   */

#include <stdio.h>
#include <stdlib.h>
#include "set.h"
#include "poly.h"
#include "qhull.h"
#include "globals.h"
#include "mem.h"

/* globals and constants */

#define MAXarraysize 8        /* max. size of freelist array */
#define BUFsize 25000         /* size of memory allocation buffer */
setT *freelist[MAXarraysize]; /* array of freelists, compiler sets to NULL */
                              /*   linked by set.e[0] */
int sizetable[MAXarraysize];  /* size of freelists, compiler sets to 0 */
int *indextable;              /* table of freelists, indexed by size */
int TABsize=0;                /* actual size of sizetable[] and freelist[] */
int MAXindex=0;               /* size of indextable (max.+1 of sizetable) */

void *buffer_begin= NULL;     /* pointer to the beginning of current buffer */

/* to avoid bus errors, memory allocation must consider alignment requirements.
   malloc() automatically takes care of alignment.   Since qhull manages
   its own memory, we need to explicitly specify alignment.
   Therefore, alignment must be explicitly specified, and that is the meaning
   of alignment_requirement -variable. sizeof(type) might be too strict.
   If gcc is available, use __alignof__(type) to get the best possible
   alignment.
   */
int alignment_requirement= sizeof(realT) - 1;     /* must be 2^n-1 */
     
/* internal functions */
  
void addsizetable(int size);
static int intcompare(const void *i, const void *j);

/* functions in alphabetical order */

/*-------------------------------------------------
-addsizetable- adds a value to sizetable if needed
*/
void addsizetable(int size) {
  int k;

  size= (size + alignment_requirement) & ~alignment_requirement;
  for(k= 0; k < TABsize; k++)
    if (sizetable[k] == size)
      return;
  if (TABsize < MAXarraysize)
    sizetable[TABsize++]= size;
  else
    fprintf(ferr, "qhull warning #7: sizetable has room only for max. %d elements\n", MAXarraysize);
} /* addsizetable */


/*-------------------------------------------------
-constructsizetable- constructs a table containing frequently used sizes
*/
void constructsizetable(void) {
  int i, k;

  if (alignment_requirement & ~alignment_requirement) {
    fprintf(ferr, "qhull error #15: alignment_requirement is not 2^n-1\n");
    errexit (NULL, NULL, NULL, NULL);
  }
  addsizetable(sizeof(vertexT));
  addsizetable(sizeof(ridgeT));
  addsizetable(sizeof(facetT));
  addsizetable(sizeof(hashtableT));
  i= sizeof(setT) + (hull_dim - 1) * SETelemsize;  /* ridge.vertices */
  addsizetable(i);
  addsizetable(normal_size);        /* normal */
  i += SETelemsize;                 /* facet.vertices, .ridges, .neighbors */
  addsizetable(i);
  qsort(sizetable, TABsize, sizeof(int), intcompare);
  MAXindex= sizetable[TABsize-1]+1;
  if (!(indextable= (int *)malloc(MAXindex * sizeof(int)))) {
    fprintf(ferr, "qhull error #14: insufficient memory\n");
    errexit(NULL, NULL, NULL, NULL);
  }
  for(k= 0; k < MAXindex; k++)
    indextable[k]= k;
  i= 0;
  for(k= 0; k < MAXindex; k++) {
    if (indextable[k] <= sizetable[i])
      indextable[k]= i;
    else
      indextable[k]= ++i;
  }
  memfree(memalloc(sizeof(ridgeT), &i), sizeof(ridgeT)); /*to initialize */
} /* constructsizetable */


/*-------------------------------------------------
-intcompare- used by qsort and bsearch to compare two integers
*/
static int intcompare(const void *i, const void *j) {
  return(*((int *)i) - *((int *)j));
} /* intcompare */


/*-------------------------------------------------
-memalloc- allocates memory for object
 uses static variables free_buffer and buffer_size
returns:
 pointer to allocated memory
 outsize= actually size allocated
*/
void *memalloc(int insize, int *outsize) {
  static void *free_buffer;    /* pointer to available memory */
  static int buffer_size;      /*  remaining bytes in free_buffer */
  setT **freelistp;
  int index;
  void *object;
  
  if (insize < MAXindex) {
    index= indextable[insize];
    *outsize= sizetable[index];
    freelistp= freelist + index;
    if (*freelistp) {
      object= (void *)(*freelistp);
      *freelistp= (*freelistp)->e[0];
      return (object);
    }else {
      if (*outsize > buffer_size) {
	object= buffer_begin;
	if (!(buffer_begin= free_buffer= malloc(BUFsize))) {
	  fprintf(ferr, "qhull error #14: insufficient memory\n");
	  buffer_begin= object;
	  errexit(NULL, NULL, NULL, NULL);
	} 
	*((void **)free_buffer)= object;
	free_buffer = buffer_begin + sizeof(realT);
	buffer_size= BUFsize - sizeof(realT); 
      }
      object= free_buffer;
      free_buffer += *outsize;
      buffer_size -= *outsize;
      return(object);
    }
  }else {
    *outsize= insize;
    if (!(object= malloc(insize))) {
      fprintf(ferr, "qhull error #14: insufficient memory\n");
      errexit(NULL, NULL, NULL, NULL);
    }
    return (object);
  }
} /* memalloc */


/*-------------------------------------------------
-memfree- frees the memory pointed by pointer
*/
void memfree(void *pointer, int size) {
  int index;
  setT *freeblock;  /* free memory blocks are all casted as setT structures */
  
  if (size < MAXindex) {
    index= indextable[size];
    freeblock= (setT *)pointer;
    freeblock->e[0]= freelist[index];
    freelist[index]= freeblock;
  }else 
    free (pointer);
} /* memfree */

