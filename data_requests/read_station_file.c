// Reads the station file filename into the list of stations station_list
//
// The station file format is as follows
/*
FDSN Station List (sorted by station)

 NETWORKSTATION LAT        LON         ELEVATIONSITE
 IU     AAE     9.0292     38.7656     2442     Addis Ababa, Ethiopia
 II     AAK     42.639     74.494      1645     Ala Archa, Kyrgyzstan
 II     ABKT    37.9304    58.1189     678      Alibek, Turkmenistan
 IU     ADK     51.8837    -176.6844   116      Adak, Aleutian Islands, Alaska
 IU     AFI     ...
*/

#include"header.h"

Station * read_station_file (char *filename, int *n_stations){

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

    sscanf(line,"%2s %5s %f %f %f %120c",station_list[i].network, station_list[i].station, &station_list[i].lat, &station_list[i].lon, &station_list[i].elevation, station_list[i].location);

    i++;

  }

  *n_stations=i;
  station_list=realloc(station_list, *n_stations*sizeof(Station));

  fclose(fp);

  return station_list;
}
