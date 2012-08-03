/* set.h -- header file for set.c that implements set
   operations needed in quickhull - application
   set operations satisfy the following properties:
    - every set is NULL terminated
    - sets may be sorted or unsorted, the caller must distinguish this
   
   Author:   Hannu Huhdanpaa    
   Internet: hannu@geom.umn.edu
    
   The Geometry Center
   
   History:

   02/06/93 - changed the FOREACH-loops (bb & hh)
   01/30/93 - setadd, setdel, setnewadd, setnewdel, setequal added (hh)
   01/25/93 - changed setmerge to use pointers (hh)
   01/21/93 - started to write (hh)
   */

/* ----------------------------------------------
-type names and included files, structures defined below
*/

typedef struct setT setT;   /* a set is a sorted array of pointers */

/* ----------------------------------------------
-constants and flags
*/
#define SETelemsize sizeof(void *) /* specifies size of set element in bytes */
#define MAX_SIZE 65535

/* ================= -structures- ===============
   see typedefs above
*/

/* ----------------------------------------------
-setT- a set of anything
*/

struct setT {
  unsigned int size;   /* MAXIMUM number of elements (except NULL) */
  void *e[1];          /* array of pointers, tail is NULL */
                       /* last slot (unless NULL) is actual size+1 */
};


/* =========== -macros- ========================= */

/*-----------------------------------------------
-FOREACHxxx- standard for loops.
        variable is NULL at end of loop
*/

#define FOREACHsetelement_(type, set, variable) variable= NULL; \
        if (set) for(\
          variable##p= (type **)&((set)->e[0]); \
	  variable= *variable##p; \
          variable##p++)
/* FOREACHsetelementreverse12 returns e[1], e[0], e[2], e[3], ... */
#define FOREACHsetelementreverse12_(type, set, variable) variable= NULL; \
        if (set) for(\
          variable##p= (type **)&((set)->e[1]); \
	  variable= *variable##p; \
          (variable##p == (type **)&((set)->e[0]))?variable##p += 2:(variable##p == (type **)&((set)->e[1]))?variable##p--:variable##p++) 
#define SETelem_(n, type, set) ((set) ? ((type *)((set)->e[n])) : NULL)
#define SETelemaddr_(n, type, set) ((set) ? ((type **)(&((set)->e[n]))) : NULL)
#define SETempty_(set) (SETelem_(0, void, set) ? False: True)
#define SETlast_(set, setsize) (((set) && !SETempty_(set))?SETelem_(SETreturnsize_(set, setsize)-1, void, set):NULL) 
#define SETreturnsize_(set, setsize) (((setsize)= (int)((set)->e[(set)->size]))?(--setsize):(set)->size)
/* SETsize- returns number of elements in set, set can't be NULL */
#define SETsize_(set, setsize) {int siz= (set)->size; \
			     if ((setsize)= (int)((set)->e[siz]))(setsize)--; \
                             else (setsize)= siz;}
#define SETdellast_(set, setsize) {\
    if((setsize)= (int)((set)->e[(set)->size])){\
      (set)->e[setsize - 2]= NULL;\
      ((int)((set)->e[(set)->size]))--;\
    }else {\
      (set)->e[(set)->size - 1]= NULL;\
      (set)->e[(set)->size]= (void*)((set)->size);\
    }}
						     

/* =========== -functions and procedures- ======= */

/********* functions in alphabetical order ********/

void setadd(setT **setp, void *elem);
void setappend(setT **setp, void *elem);
void setappend2ndlast(setT **setp, void *elem);
setT *setcopy(setT *set);
void setdel(setT **setp, void *elem);
int setequal(setT *setA, setT *setB);
void setfree(setT **set);
int setin(setT *set, void *element, int upper);
void setlarger(setT **setp);
setT *setnew(int *setsize);
setT *setnewdel(setT *set, void *elem);
setT *setnewprepend(setT *set, void *newelem);
void setprepend(setT **setp, void *elem);
void setreplace(setT **setp, void *oldelem, void *newelem);
int setsize(setT *size);




