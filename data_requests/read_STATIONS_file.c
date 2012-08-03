// Reads the station file filename into the list of stations station_list
//
// The station file format is as follows
/*
110
  ADK  IU   51.8837  -176.6844   116.0   0.0
  AFI  IU  -13.9093  -171.7773   706.0   0.0
 ANMO  IU   34.9462  -106.4567  1840.0   0.0
 BILL  IU   68.0651   166.4524   299.0   0.0

*/

#include"header.h"

Station * read_STATIONS_file (char *filename, int *n_stations){

  FILE *fp;
  char line[LINE_LEN];
  int i, list_length;
  Station *station_list;

  fp=fopen(filename, "r");

  fscanf(fp,"%d",&list_length);
  // Allocate memory for stations list
  station_list=malloc(list_length*sizeof(Station));

  // Read each subsequent line up to the end of the file and
  // add to station list
  i=0;
  for (i=0;i<list_length;i++){
      fscanf(fp,"%5c %2c %f %f %f %f",station_list[i].station, station_list[i].network, &station_list[i].lat, &station_list[i].lon, &station_list[i].elevation, &station_list[i].depth);
      fprintf(stdout,"%5s %2s %f %f %f %f\n",station_list[i].station, station_list[i].network, station_list[i].lat, station_list[i].lon, station_list[i].elevation, station_list[i].depth);
  }

  *n_stations=list_length;

  fclose(fp);

  return station_list;
}
