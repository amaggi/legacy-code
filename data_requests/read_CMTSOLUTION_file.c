// Reads the event file filename into the list of events event_list
//
// The CMTSOLUTION file format is 
/*
PDE 1994  1 17 12 30 55.30  34.2100 -118.5400  18.0 6.4 6.8 SOUTHERN CALIFORNIA                     
event name:     011794B
time shift:      0.0000
half duration:   0.0000
latitude:       34.4400
longitude:    -118.6400
depth:          16.8000
Mrr:       1.081000e+26
Mtt:      -9.430000e+25
Mpp:      -1.380000e+25
Mrt:       5.400000e+24
Mrp:      -3.990000e+25
Mtp:       4.360000e+25
*/


#include"header.h"

Event *read_CMTSOLUTION_file (char *filename, int *n_events){

  FILE *fp;
  char line[LINE_LEN];
  int i, list_length;
  Event *event_list;

  // Allocate event list to be 1 event long
  list_length=1;
  event_list=malloc(list_length*sizeof(Event));

  fp=fopen(filename, "r");

  // Read the first line (header of file)
  fgets(line, LINE_LEN, fp);
  sscanf(line,"%*s %d %d %d %d %d %f %f %f %f %f",&event_list[i].year, &event_list[i].month, &event_list[i].day, &event_list[i].hour, &event_list[i].min, &event_list[i].sec, &event_list[i].lat, &event_list[i].lon, &event_list[i].depth, &event_list[i].mag); 

//  fprintf(stderr,"%d %d %d %d %d %f %f %f %f %f\n",event_list[i].year, event_list[i].month, event_list[i].day, event_list[i].hour, event_list[i].min, event_list[i].sec, event_list[i].lat, event_list[i].lon, event_list[i].depth, event_list[i].mag); 
  // get the second line (with the identifier)
  fgets(line, LINE_LEN, fp);
  sscanf(line,"event name: %s",event_list[i].identifier);

  *n_events=1;

  fclose(fp);

  return event_list;
}
