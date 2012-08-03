#include "header.h"

/* read and interpret commant line, and opens input/output files
 * returns 1 on success, 0 on failure */

int parser (int argc, char *argv[]){

  int i;
  int year, month, day, hour, min;

  if(argc == 1) { /* no arguments */
    print_help();
    return 0;
  }

  /* initialize all flags to not-set */
  cat_type=UNSET;
  cmt_hypo=UNSET;
  mag_search=UNSET;
  dep_search=UNSET;
  circ_search=UNSET;
  seg_search=UNSET;
  rec_search=UNSET;
  date_search=UNSET;
  output_format=OUT_TXT;
  header_format=SHORT;

  strcpy(basename,"search");

  /* step through the arguments */

  i=1;
  while(i < argc){
    /* help */
    if(strcmp(argv[i],"-help")==0){
      print_help();
      return 0;
    }
    /* catalogue type */
    else if(strcmp(argv[i],"-cat")==0){
      /* check if there are enough arguments left */
      if(argc < i+3){
        print_help();
        return 0;
      }
      i++;
      /* check for catalogue type */
      if(strcmp(argv[i],"cmt")==0){
       cat_type=CMT;
       i++;
      } 
      else if(strcmp(argv[i],"pde")==0){
       cat_type=PDE;
       i++;
      } 
      else if(strcmp(argv[i],"ehb")==0){
       cat_type=EHB;
       i++;
      } 
      else if(strcmp(argv[i],"txt")==0){
       cat_type=TXT;
       i++;
      } 
      else if(strcmp(argv[i],"jdy")==0){
       cat_type=JDY;
       i++;
      } 
      else{ /* unknown type - return 0 */
        fprintf(stderr,"Unknown catalogue type - exiting\n");
        return 0;
      }
      /* catalogue type set correctly - read filename */
      strcpy(cat_filename,argv[i]);
    }
    /* use cmt hypocenter */
    else if(strcmp(argv[i],"-cen")==0){
      cmt_hypo=SET;
    }
    /* moment search */
    else if(strcmp(argv[i],"-mom")==0){
      if(argc < i+3){
        print_help();
        return 0;
      }
      mag_search=MOM;
      i++;
      minmo=(float)atof(argv[i]);
      i++;
      maxmo=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-mag")==0){
      if(argc < i+3){
        print_help();
        return 0;
      }
      mag_search=MAG;
      i++;
      minmag=(float)atof(argv[i]);
      i++;
      maxmag=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-ms")==0){
      if(argc < i+3){
        print_help();
        return 0;
      }
      mag_search=MS;
      i++;
      minmag=(float)atof(argv[i]);
      i++;
      maxmag=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-mb")==0){
      if(argc < i+3){
        print_help();
        return 0;
      }
      mag_search=MB;
      i++;
      minmag=(float)atof(argv[i]);
      i++;
      maxmag=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-dep")==0){
      if(argc < i+3){
        print_help();
        return 0;
      }
      dep_search=SET;
      i++;
      mindep=(float)atof(argv[i]);
      i++;
      maxdep=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-cir")==0){
      if(argc < i+5){
        print_help();
        return 0;
      }
      circ_search=SET;
      i++;
      cenlat=(float)atof(argv[i]);
      i++;
      cenlon=(float)atof(argv[i]);
      i++;
      rmin=(float)atof(argv[i]);
      i++;
      rmax=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-rec")==0){
      if(argc < i+5){
        print_help();
        return 0;
      }
      rec_search=SET;
      i++;
      minlat=(float)atof(argv[i]);
      i++;
      maxlat=(float)atof(argv[i]);
      i++;
      minlon=(float)atof(argv[i]);
      i++;
      maxlon=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-seg")==0){
      if(argc < i+7){
        print_help();
        return 0;
      }
      seg_search=SET;
      i++;
      cenlat=(float)atof(argv[i]);
      i++;
      cenlon=(float)atof(argv[i]);
      i++;
      rmin=(float)atof(argv[i]);
      i++;
      rmax=(float)atof(argv[i]);
      i++;
      bazmin=(float)atof(argv[i]);
      i++;
      bazmax=(float)atof(argv[i]);
    }
    else if(strcmp(argv[i],"-date")==0){
      if(argc < i+5){
        print_help();
        return 0;
      }
      date_search=SET;
      start_time.tm_isdst=-1; /* avoid daylight saving issues */
      end_time.tm_isdst=-1;
      i++;
      sscanf(argv[i],"%4d%2d%2d",&year, &month, &day);
      start_time.tm_year = year - 1900;
      start_time.tm_mon = month -1;
      start_time.tm_mday = day;
      i++;
      sscanf(argv[i],"%2d%2d",&hour,&min);
      start_time.tm_hour = hour;
      start_time.tm_min = min;
      i++;
      sscanf(argv[i],"%4d%2d%2d",&year,&month,&day);
      end_time.tm_year = year - 1900;
      end_time.tm_mon = month -1;
      end_time.tm_mday = day;
      i++;
      sscanf(argv[i],"%2d%2d",&hour,&min);
      end_time.tm_hour = hour;
      end_time.tm_min = min;

      mintime=mktime(&start_time);
      maxtime=mktime(&end_time);

    }
    else if(strcmp(argv[i],"-noh")==0){
      header_format=NONE;
    }
    else if(strcmp(argv[i],"-lhd")==0){
      header_format=LONG;
    }
    else if(strcmp(argv[i],"-jdy")==0){
      output_format=OUT_JDY;
    }
    else if(strcmp(argv[i],"-reg")==0){
      output_format=OUT_REG;
    }
    else if(strcmp(argv[i],"-out")==0){
      if(argc < i+2){
        print_help();
        return 0;
      }
      i++;
      strcpy(basename,argv[i]);
    }
    else{
      print_help();
      return 0;
    }
    i++;
  }

  if(cat_type == UNSET) {
    fprintf(stderr,"Need to specify a catalogue!\n");
    return 0;
  }

  return 1;
}


void print_parsed_info(){

  char timestring[80];

  printf("Parsed information:\n");
  printf("\n");
  printf("Catalogue type: ");
  switch (cat_type) {
    case CMT:
      printf("CMT\n");
      break;
    case PDE:
      printf("PDE\n");
      break;
    case EHB:
      printf("EHB\n");
      break;
    case TXT:
      printf("TXT\n");
      break;
    case JDY:
      printf("JDY\n");
      break;
    default :
      printf("Error!\n");
      break;
  }
  printf("Catalogue filenamename: %s\n",cat_filename);
  printf("Output basename: %s\n",basename);
  if(cmt_hypo == SET) printf("Use CMT hypocenters\n");
  if(dep_search == SET) printf("Depth search: %.2f %.2f\n",mindep, maxdep);
  if(circ_search == SET) printf("Circle search: %.2f %.2f %.2f %.2f\n",cenlat, cenlon, rmin, rmax);
  if(seg_search == SET) printf("Segment search: %.2f %.2f %.2f %.2f %.2f %.2f\n",cenlat, cenlon, rmin, rmax, bazmin, bazmax);
  if(rec_search == SET) printf("Rectangle search: %.2f %.2f %.2f %.2f\n",minlat, maxlat, minlon, maxlon);
  if(date_search == SET) {
    printf("Date search: ");
    strftime(timestring,80,"%Y-%m-%d %H-%M",&start_time);
    printf("%s to ",timestring);
    strftime(timestring,80,"%Y-%m-%d %H-%M",&end_time);
    printf("%s\n",timestring);
  }
  if(mag_search != UNSET){
    printf("Size search - ");
    switch (mag_search){
      case MOM:
        printf("moment : %.2f %.2f\n", minmo, maxmo);
        break;
      case MAG:
        printf("magnitude : %.2f %.2f\n", minmag, maxmag);
        break;
      case MS:
        printf("ms : %.2f %.2f\n", minmag, maxmag);
        break;
      case MB:
        printf("mb : %.2f %.2f\n", minmag, maxmag);
        break;
      default:
        printf("Error!\n");
        break;
      break;
    }
  }
  printf("Output format : ");
  switch(output_format){
    case OUT_TXT:
      printf("old style .txt\n");
      break;
    case OUT_REG:
      printf("old style .txt with region \n");
      break;
    case OUT_JDY:
      printf("new style .txt (jdy)\n");
      break;
    default:
      printf("Error!\n");
      break;
  }
  printf("Header format : ");
  switch(header_format){
    case SHORT:
      printf("short\n");
      break;
    case LONG:
      printf("lhd\n");
      break;
    case NONE:
      printf("no header\n");
      break;
    default:
      printf("Error!\n");
      break;
  }
  printf("\n");
}
