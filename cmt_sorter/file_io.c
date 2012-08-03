#include"header.h"


struct PDE_Record pde;
struct CMT_Record cmt;
struct EHB_Record ehb;
struct JDY_Record jdy;

void print_headers(int argc, char *argv[]){

  int i;
  
  if(header_format == NONE)
    return;

  if(header_format == LONG){
    for(i=0;i<argc;i++)
      fprintf(fpout,"%s ", argv[i]);
    fprintf(fpout,"\n");
  }

  if(output_format == OUT_JDY) 
    fprintf(fpout,"          Origin time           Lat      Lon     Dp  Mag\n");

  else { 
    fprintf(fpout,"                             PDE                      CMT                      NP1         NP2\n");

    if(output_format == OUT_REG)
      fprintf(fpout,"          Origin time        Lat    Long   Dep   mb  Ms     Lat    Long    Dep Scalar Mmt   Str Dip  Rk  Str Dip  Rk  Region                 name\n");

    if(output_format == OUT_TXT)
      fprintf(fpout,"  DATE  TIME        Lat    Long   Dep   mb  Ms     Lat    Long    Dep Scalar Mmt   Str Dip  Rk  Str Dip  Rk\n");

  }

}

int read_event(){


  /* null all cmt specific stuff */
  event.cmt_lat=0;
  event.cmt_lon=0;
  event.cmt_dep=0;
  event.mom_base=0;
  event.mom_exp=0;
  event.strike1=0;
  event.dip1=0;
  event.rake1=0;
  event.strike2=0;
  event.dip2=0;
  event.rake2=0;

  /* now read according to type */
  switch(cat_type){
    case CMT:
      read_cmt_record();
      assign_cmt_event();
      break;
    case PDE:
      read_pde_record();
      assign_pde_event();
      break;
    case EHB:
      read_ehb_record();
      assign_ehb_event();
      break;
    case TXT:
      break;
    case JDY:
      read_jdy_record();
      assign_jdy_event();
      break;
    default:
      return 0;
      break;
  }

  return 1;

}

void write_event(){

  char evname[14];
  char datestring[80];
  float mag;

  /*
  if(mag_search == MB)
    mag=event.mb;
  else
    mag=event.ms;
  */
  mag=event.mb;
  
  switch(output_format){

    case OUT_TXT:
      strftime(datestring,80,"%Y%m%d %H%M%S", &event.origin_time);
      fprintf(fpout,"%s.%d %6.2f %7.2f %6.2f %3.1f %3.1f %6.2f %7.2f %6.2f %5.2f e%3d %4d %3d %4d %4d %3d %4d\n",datestring, event.msec/100, event.lat, event.lon, event.dep, event.mb, event.ms, event.cmt_lat, event.cmt_lon, event.cmt_dep, event.mom_base, event.mom_exp, event.strike1, event.dip1, event.rake1, event.strike2, event.dip2, event.rake2);
      break;

    case OUT_JDY:
      set_tm_jday(&event.origin_time);
      strftime(datestring,80,"%Y (%j) %m %d %H %M %S", &event.origin_time);
      strftime(evname,14,"%Y%j%H%M%S", &event.origin_time);
      fprintf(fpout,"%s.%d  %7.3f  %8.3f %4d  %3.1f %s\n",datestring,event.msec/100,event.lat, event.lon, (int)event.dep,mag,evname);
      break;

    case OUT_REG:
      set_tm_jday(&event.origin_time);
      strftime(datestring,80,"%Y (%j) %m %d %H %M %S", &event.origin_time);
      strftime(evname,14,"%Y%j%H%M%S", &event.origin_time);
      fprintf(fpout,"%s.%d %6.2f %7.2f %6.2f %3.1f %3.1f %6.2f %7.2f %6.2f %5.2f e%3d %4d %3d %4d %4d %3d %4d %s %s\n",datestring,event.msec/100,event.lat, event.lon, event.dep, event.mb, event.ms, event.cmt_lat, event.cmt_lon, event.cmt_dep, event.mom_base, event.mom_exp, event.strike1, event.dip1, event.rake1, event.strike2, event.dip2, event.rake2, event.region, evname);
      break;
  }
}

void read_ehb_record(){

  char dummy;
  
  fread(&ehb, sizeof(struct EHB_Record),1,fpin);
  if(!feof(fpin))
    fread(&dummy, sizeof(char),1,fpin); 

}

void assign_ehb_event(){

  int year, month;

  year=fixed_len_atoi(ehb.yr,2);

  if (year < 50)
    year = 2000 + year;
  else
    year = 1900 + year;

  month=fixed_len_atoi(ehb.mon,2);

  event.origin_time.tm_year = year - 1900;
  event.origin_time.tm_mon = month - 1;
  event.origin_time.tm_mday = fixed_len_atoi(ehb.day,2);
  event.origin_time.tm_hour = fixed_len_atoi(ehb.hr,2);
  event.origin_time.tm_min = fixed_len_atoi(ehb.min,2);
  event.origin_time.tm_sec = (int)fixed_len_atof(ehb.sec,6);
  event.msec=(int)((fixed_len_atof(ehb.sec,6)-event.origin_time.tm_sec)*1000);

  event.lat=(float)fixed_len_atof(ehb.lat,8);
  event.lon=(float)fixed_len_atof(ehb.lon,8);
  event.dep=(float)fixed_len_atof(ehb.dep,6);
  event.mb=(float)fixed_len_atof(ehb.dmag,4);
  event.ms=(float)fixed_len_atof(ehb.imag,2);

}


void read_cmt_record(){

  char dummy;
  
  fread(&cmt, sizeof(struct CMT_Record),1,fpin);
  if(!feof(fpin))
    fread(&dummy, sizeof(char),1,fpin); 

}

void assign_cmt_event(){

  int year, month;

  year=fixed_len_atoi(cmt.two_digit_year,2);

  if (year < 50)
    year = 2000 + year;
  else
    year = 1900 + year;

  month=fixed_len_atoi(cmt.month,2);

  event.origin_time.tm_year = year - 1900;
  event.origin_time.tm_mon = month - 1;
  event.origin_time.tm_mday = fixed_len_atoi(cmt.day,2);
  event.origin_time.tm_hour = fixed_len_atoi(cmt.hour,2);
  event.origin_time.tm_min = fixed_len_atoi(cmt.min,2);
  event.origin_time.tm_sec = (int)fixed_len_atof(cmt.second,4);
  event.msec=(int)((fixed_len_atof(cmt.second,4)-event.origin_time.tm_sec)*1000);


  event.lat=(float)fixed_len_atof(cmt.lat,6);
  event.lon=(float)fixed_len_atof(cmt.lon,7);
  event.dep=(float)fixed_len_atof(cmt.depth,5);
  event.mb=(float)fixed_len_atof(cmt.mb,3);
  event.ms=(float)fixed_len_atof(cmt.ms,3);

  event.cmt_lat=(float)fixed_len_atof(cmt.cmt_lat,6);
  event.cmt_lon=(float)fixed_len_atof(cmt.cmt_lon,7);
  event.cmt_dep=(float)fixed_len_atof(cmt.cmt_dep,5);
  event.mom_base=(float)fixed_len_atof(cmt.cmt_base,4);
  event.mom_exp=fixed_len_atoi(cmt.cmt_mom,2);

  event.strike1=fixed_len_atoi(cmt.strike1,3);
  event.dip1=fixed_len_atoi(cmt.dip1,2);
  event.rake1=fixed_len_atoi(cmt.rake1,4);
  event.strike2=fixed_len_atoi(cmt.strike2,3);
  event.dip2=fixed_len_atoi(cmt.dip2,2);
  event.rake2=fixed_len_atoi(cmt.rake1,4);

  strncpy(event.region,cmt.region,24);

}

void read_pde_record(){

  char dummy;

  fread(&pde, sizeof(struct PDE_Record),1,fpin);
  
  /* read the end of line carriage return - beware, may not be one at end of
 * last line in file */

  if(!feof(fpin)) 
    fread(&dummy, sizeof(char),1,fpin); 

}

void assign_pde_event(){

  int year, month;

  year=fixed_len_atoi(pde.Year,4);
  month=fixed_len_atoi(pde.Month,2);
  event.origin_time.tm_year = year - 1900;
  event.origin_time.tm_mon = month - 1;
  event.origin_time.tm_mday = fixed_len_atoi(pde.Day,2);
  event.origin_time.tm_hour = fixed_len_atoi(pde.Hour,2);
  event.origin_time.tm_min = fixed_len_atoi(pde.Min,2);
  event.origin_time.tm_sec = (int)fixed_len_atof(pde.Sec,5);
  event.msec=(int)((fixed_len_atof(pde.Sec,5)-event.origin_time.tm_sec)*1000);

  event.lat=(float)fixed_len_atof(pde.Lat,7);
  event.lon=(float)fixed_len_atof(pde.Lon,8);
  event.dep=(float)fixed_len_atof(pde.Depth,3);
  event.mb=(float)fixed_len_atof(pde.mb,3);
  event.ms=(float)fixed_len_atof(pde.Ms,3);
  
}

void read_jdy_record(){

  char dummy;

  fread(&jdy, sizeof(struct JDY_Record),1,fpin);
  
  /* read the end of line carriage return - beware, may not be one at end of
 * last line in file */

  if(!feof(fpin)) 
    fread(&dummy, sizeof(char),1,fpin); 

}

void assign_jdy_event(){

  int year, month;

  year=fixed_len_atoi(jdy.year,4);
  month=fixed_len_atoi(jdy.month,2);
  event.origin_time.tm_year = year - 1900;
  event.origin_time.tm_mon = month - 1;
  event.origin_time.tm_mday = fixed_len_atoi(jdy.day,2);
  event.origin_time.tm_hour = fixed_len_atoi(jdy.hour,2);
  event.origin_time.tm_min = fixed_len_atoi(jdy.min,2);
  event.origin_time.tm_sec = (int)fixed_len_atof(jdy.sec,4);

  event.lat=(float)fixed_len_atof(jdy.lat,7);
  event.lon=(float)fixed_len_atof(jdy.lon,8);
  event.dep=(float)fixed_len_atof(jdy.depth,3);
  event.mb=(float)fixed_len_atof(jdy.mag,3);
  event.ms=0;
  
}


