#include"header.h"

int check_criteria(){

  int satisfy=1;
  float lat, lon, dep;
  float elat, elon, slat, slon, azm, bzm, ddg, dkm;
  float delta1, delta2;
  time_t origin;

  /* use cmt values if required */
  if (cmt_hypo == SET){
    dep=event.cmt_dep;
    lat=event.cmt_lat;
    lon=event.cmt_lon;
  }
  else {
    dep=event.dep;
    lat=event.lat;
    lon=event.lon;
  }

  /* if ms = 0 then set ms = mb for purpose of searching */
  if (event.ms == 0)
   event.ms = event.mb;

  /* depth search */
  if (dep_search == SET) 
    if(dep < mindep || dep > maxdep)
      satisfy *= 0;
  
  /* rectangle search */
  if (rec_search == SET) 
    if(lat < minlat || lat > maxlat || lon < minlon || lon > maxlon)
      satisfy *= 0;
 
  /* magnitude search */
  switch(mag_search){
    case UNSET:
      /* do nothing */
      break;
    case MOM:
      if (event.mom_exp < minmo || event.mom_exp > maxmo)
        satisfy *= 0;
      break;
    case MAG:
      /* fall into case MS */
    case MS:
      if (event.ms < minmag || event.ms > maxmag)
        satisfy *= 0;
      break;
    case MB:
      if (event.mb < minmag || event.mb > maxmag)
        satisfy *= 0;
      break;
  }

  /* date search */
  if(date_search == SET){
    origin=mktime(&event.origin_time);
    if(difftime(origin, maxtime) > 0 || difftime(mintime, origin) > 0)
        satisfy *= 0;
  }

  /* circle search */
  if(circ_search == SET){
    slat=cenlat;
    slon=cenlon;
    elat=lat;
    elon=lon;
    distaz_(&elat, &elon, &slat, &slon, &azm, &bzm, &ddg, &dkm);
    if(ddg < rmin || ddg > rmax)
        satisfy *= 0;
  }

   /* segment search */
// Fix this at a later date!
  if(seg_search == SET){
    slat=cenlat;
    slon=cenlon;
    elat=lat;
    elon=lon;
    distaz_(&elat, &elon, &slat, &slon, &azm, &bzm, &ddg, &dkm);

//    delta1=360-bazmin;
//    delta2=sum_angle(bazmax, delta1);
//    bzm=sum_angle(bzm,delta1);

    if(ddg < rmin || ddg > rmax || bzm < bazmin || bzm > bazmax)
        satisfy *= 0;
  }
  

  return satisfy;

}
