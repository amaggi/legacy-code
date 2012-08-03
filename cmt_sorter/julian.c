#include<stdio.h>
#include<time.h>

extern void julian_(int *year, int *month, int *day, int *jday);

void set_tm_jday(struct tm *tm){
 
  int year, month, day, jday;
  
  year=tm->tm_year+1900;
  month=tm->tm_mon+1;
  day=tm->tm_mday;

  julian_(&year, &month, &day, &jday);

  tm->tm_yday=jday-1;

}

