/*
* $Id: maxima.c,v 1.1 2005/03/24 03:52:41 alessia Exp $
*
* Given an array of floating point numbers, returns an array of 
* integers containing the indexes of the maxima of the array
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "amutils.h"

#define UP 1
#define DOWN -1
#define FLAT 0

#define START_MAX 10

#define EPS 1.0E-9

int direction(float a, float b, float w_level);

void maxima(float *x, int nx, int *maxima , int *nindex, float w_level)
{
   int i,k;
   int prev_dir, next_dir;
   
   k=0;

/* set up the starting previous direction */
   prev_dir=direction(x[0],x[1],w_level);

/* iterate through the points in the array */
   for (i=1;i<nx-1;i++) {
/*   the next direction */
     next_dir=direction(x[i],x[i+1],w_level);
     if (next_dir==DOWN && prev_dir == UP) {
       maxima[k]=i;
       k++;
     }
/*   set up prev_dir for next iternation */
     prev_dir = next_dir;
   }

   (*nindex)=k;
}

int direction(float a, float b, float w_level) {
  if ( a < w_level || b < w_level ) return FLAT;
  if ( fabs(a-b) < EPS) return FLAT;
  if (b > a) return UP; 
  else return DOWN;
}

/*given an array, the index of a local maximum, and a water level, give the index of the closest points in the array that dip below the water level*/
void bracket_maximum(float *x, int nx, int imax , float w_level, int *i1, int *i2) {

  int i;
  
  (*i1)=0;
  (*i2)=nx-1;

  for(i=imax;i>=0;i--){
    if(x[i]<w_level){
      (*i1) = i ; 
      break;
    }
  }
  for(i=imax;i<nx;i++){
    if(x[i]<w_level){
      (*i2) = i ; 
      break;
    }
  }

}

