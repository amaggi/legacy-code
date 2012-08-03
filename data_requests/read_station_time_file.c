// Reads the station file filename into the list of stations station_list
//
// The station file format is as follows
/*
IRIS Webquery: Station Information

NETWORK STATION STARTTIME       ENDTIME LAT     LON     ELEVATION
AK      MCK     2000-06-08 00:00:00.0   2015-12-31 23:59:59.0   +63.732300      -148.934900     618
AL      A011    1966-01-01 00:00:00.0   1979-01-01 00:00:00.0   +65.233330      -147.743330     640
AL      A021    1966-01-01 00:00:00.0   1979-01-01 00:00:00.0   +65.373610      -147.401110     564

*/

#include"header.h"

Station * read_station_time_file (char *filename, int *n_stations){

  FILE *fp;
  char line[LINE_LEN];
  int i, list_length;
  Station *station_list;

  fp=fopen(filename, "r");

  // Allocate station list to be N_STATIONS long
  list_length = N_STATIONS;
  station_list=malloc(list_length*sizeof(Station));

  // Read the first 3 lines (header of file)
  fgets(line, LINE_LEN, fp);
  fgets(line, LINE_LEN, fp);
  fgets(line, LINE_LEN, fp);

  // Read each subsequent line up to the end of the file and
  // add to station list
  i=0;
  while(fgets(line, LINE_LEN, fp)!=NULL){
 
    if(i == list_length){
      // have reached current limit of station_list
      list_length*=2;
      station_list=realloc(station_list, list_length*sizeof(Station));
    }

    sscanf(line,"%2s %5s %*s %*s %*s %*s %f %f %f",station_list[i].network, station_list[i].station, &station_list[i].lat, &station_list[i].lon, &station_list[i].elevation);

    i++;

  }

  *n_stations=i;
  station_list=realloc(station_list, *n_stations*sizeof(Station));

  fclose(fp);

  return station_list;
}
