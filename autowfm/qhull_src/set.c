/* set.c -- implements set manipulations needed for quickhull 
   see set.h for further details

   Author: Hannu Huhdanpaa
   Email : hannu@geom.umn.edu

   The Geometry Center
*/

#include <stdio.h>
#include <memory.h>
#include "mem.h"
#include "set.h"

#define FOREACHelem_(set) FOREACHsetelement_(void, set, elem)
#define addto_(p, list)   {(p)->next= (list); (list)= (p);}

/*----------------------------------------
-setadd- adds an element to a sorted set
  setp may be NULL, newlem must not be NULL
*/
void setadd(setT **setp, void *newelem) {
  int size;
  setT *newset;
  void **elemOp, **elemNp;

  if (!(*setp)) {
    size= 1;
    *setp= setnew(&size);
  }
  SETsize_(*setp, size);
  if (setin(*setp, newelem, size-1))
    return;
  if (size == (*setp)->size)
    setlarger(setp);
  size++;
  (*setp)->e[(*setp)->size]= (void *)(size + 1); /* 1 relative */
  elemNp= &((*setp)->e[size]);
  *elemNp--= NULL;
  elemOp= &((*setp)->e[size-2]);
  if (size == 1 || (*setp)->e[0] > newelem) {
    void **elem0p= &(*setp)->e[0];
    while (elemOp >= elem0p)
      *elemNp--= *elemOp--;
  }else {
    while(*elemOp > newelem)
      *elemNp--= *elemOp--;
  }
  *elemNp= newelem;
} /* setadd */


/*----------------------------------------
-setappend- appends an element into an unsorted set
newelem is not allowed to be NULL
*/
void setappend(setT **setp, void *newelem) {
  int size;
  setT *newset;
  void **elemOp, **elemNp;
  int *sizep;
  void **endp;
  
  if (!(*setp)) {
    size= 1;
    *setp= setnew(&size);
  }
  sizep= (int *)(&((*setp)->e[(*setp)->size]));
  if (!*sizep) {
    setlarger(setp);
    sizep= (int *)(&((*setp)->e[(*setp)->size]));
  }
  *(endp= &((*setp)->e[(*sizep)++ - 1]))= newelem;
  *(++endp)= NULL;
} /* setappend */


/*----------------------------------------
-setappend2ndlast- appends an element next to the last element into an 
 unsorted set
newelem is not allowed to be NULL
*/
void setappend2ndlast(setT **setp, void *newelem) {
  int size;
  setT *newset;
  void **elemOp, **elemNp;
  
  if (!(*setp)) {
    size= 1;
    *setp= setnew(&size);
  }
  SETsize_(*setp, size);
  if (size == (*setp)->size)
    setlarger(setp);
  (*setp)->e[size]= (*setp)->e[size - 1];
  (*setp)->e[size-1]= newelem;
  (*setp)->e[(*setp)->size]= (void *)(((int)(*setp)->e[(*setp)->size]) + 1);
  (*setp)->e[size+1]= NULL;
} /* setappend2ndlast */


/*----------------------------------------
-setcopy- copies a sorted or unsorted set into another
returns:
  new set has the same maximum size as the old set
*/
setT *setcopy(setT *set) {
  setT *newset;
  int size;
  
  size= set->size;
  newset= setnew(&size);
  memcpy(&(newset->e[0]), &(set->e[0]), SETelemsize * (set->size + 1));
  return (newset);
} /* setcopy */


/*----------------------------------------
-setdel- deletes newelem from sorted or unsorted set
newelm is not allowed to be NULL
*/
void setdel(setT **setp, void *newelem) {
  void **elemAp, **elemBp;
  int *sizep;

  if (!*setp)
    return;
  elemAp= SETelemaddr_(0, void, *setp);
  while(*elemAp != newelem && *elemAp)
    elemAp++;
  if (*elemAp) {
    elemBp= elemAp+1;
    while(*elemAp++= *elemBp++)
      ;
    if (*(sizep= (int *)(&((*setp)->e[(*setp)->size]))) == 1)
      *sizep= 0;
    else if (!((*sizep)--)) /*  full set */
      *sizep= (*setp)->size;  /* actual size + 1 */
  }
} /* setdel */


/*----------------------------------------
-setequal- tests whether two sorted sets are equal
*/
int setequal(setT *setA, setT *setB) {
  void **elemAp, **elemBp;
  
  elemAp= SETelemaddr_(0, void, setA);
  elemBp= SETelemaddr_(0, void, setB);
  while(*elemAp) {
    if (*elemAp++ != *elemBp++)
      return 0;
  }
  return 1;
} /* setequal */


/*----------------------------------------
-setfree- frees the space occupied by a sorted or unsorted set
*/
void setfree(setT **setp) {
  if (*setp) {
    memfree(*setp, sizeof(setT) + ((*setp)->size)*SETelemsize); 
    *setp= NULL;
  }
} /* setfree */


/*----------------------------------------
-setin- returns 1 if setelem is in a sorted set, 0 otherwise.  Internal only.
    index of last element, or -1
*/
int setin(setT *set, void *setelem, int upper) {
  int low, high, mid;
  void *elem, **elemp;

  if (upper > 5) {
    low= 0;
    high= upper;
    while (low <= high) {
      mid= (low + high)/2;
      elem= SETelem_(mid, void, set);
      if (setelem < elem)
	high= mid - 1;
      else if (setelem > elem)
	low= mid + 1;
      else
	return 1;
    }
    return 0; 
  }else {
    FOREACHelem_(set)
      if (elem == setelem)
	return 1;
    return 0;
  }
} /* setin */


/*----------------------------------------
-setlarger- returns a larger set that contains elements of *setp
*/
void setlarger(setT **setp) {
  int newsize, currentsize= (*setp)->size + 1;
  struct setT *newset;
  
  newsize= 2 * currentsize;
  newset= setnew(&newsize);
  memcpy(&(newset->e[0]), &((*setp)->e[0]), currentsize * SETelemsize);
  newset->e[newset->size]= (void *)currentsize;
  setfree(setp);
  *setp= newset;
} /* setlarger */


/*----------------------------------------
-setnew- creates and allocates space for a set
    setsize means the number of elements (NOT including the NULL terminator)
*/
setT *setnew(int *setsize) {
  setT *set;
  int sizewanted, sizereceived;
  
  sizewanted= sizeof(setT) + *setsize*SETelemsize;
  set= (setT *)memalloc(sizewanted, &sizereceived);
  *setsize += (sizereceived - sizewanted)/SETelemsize;
  set->size= *setsize;
  set->e[*setsize]= (void *)1;
  set->e[0]= NULL;             /* may overwrite end */
  return (set);
} /* setnew */


/*----------------------------------------
-setnewdel- creates a new sorted or unsorted set not containing newelem
    newelem is not allowed to be NULL
*/
setT *setnewdel(setT *set, void *newelem) {
  setT *newset;
  void **elemOp, **elemNp;
  int size;
  
  if (!set)
    return NULL;
  size= set->size;
  newset= setnew(&size);
  elemNp= SETelemaddr_(0, void, newset);
  elemOp= SETelemaddr_(0, void, set);
  while(*elemOp && *elemOp != newelem)
    *elemNp++= *elemOp++;
  if (*elemOp) {
    elemOp++;
    SETsize_(set, size);
    newset->e[newset->size]= (void *)size;  /* size-1 one relative */
  }else
    newset->e[set->size]= set->e[set->size];
  while(*elemNp++= *elemOp++)
    ;
  return(newset);
} /* setnewdel */


/*----------------------------------------
-setnewprepend- creates a new set by prepending element into it
*/
setT *setnewprepend(setT *set, void *newelem) {
  setT *newset;
  void **elemOp, **elemNp;
  int size;
  
  if (!set) {
    size= 1;
    set= setnew(&size);
  }
  SETsize_(set, size);
  size++;
  newset= setnew(&size);
  elemNp= SETelemaddr_(0, void, newset);
  elemOp= SETelemaddr_(0, void, set);
  *elemNp++= newelem;
  while(*elemNp++= *elemOp++)
    ;
  return(newset);
} /* setnewprepend */


/*----------------------------------------
-setprepend- prepends an element into a set
newelem is not allowed to be NULL
*/
void setprepend(setT **setp, void *newelem) {
  setT *newset;
  void **elemOp, **elemNp;
  int size;

  if (!(*setp)) {
    size= 1;
    *setp= setnew(&size);
  }
  SETsize_(*setp, size);
  if (size == (*setp)->size)
    setlarger(setp);
  (*setp)->e[(*setp)->size]= (void *)(((int)((*setp)->e[(*setp)->size])) + 1);
  elemNp= SETelemaddr_(size+1, void, *setp);
  elemOp= SETelemaddr_(size, void, *setp);
  for(size++; size; size--) 
    *elemNp--= *elemOp--;      
  *elemNp= newelem;
} /* setprepend */


/*----------------------------------------
-setreplace- deletes oldelem from sorted or unsorted set and inserts newelem
into its place
*/
void setreplace(setT **setp, void *oldelem, void *newelem) {
  void **elemAp, **elemBp;
  
  if (!*setp)
    return;
  elemAp= SETelemaddr_(0, void, *setp);
  while(*elemAp != oldelem && *elemAp)
    elemAp++;
  if (*elemAp)
    *elemAp= newelem;
} /* setreplace */


/*----------------------------------------
-setsize- returns the size of a set
*/
int setsize(setT *set) {
  int size;
  
  if (!set)
    return (0);
  if (size= (int)(set->e[set->size]))
    return (size-1);
  else
    return (set->size);
} /* setsize */
