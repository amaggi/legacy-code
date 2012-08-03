// Reads the event file filename into the list of events event_list
//
// The event file format is as follows
/*
           Origin time           Lat      Lon     Dp  Mag
1992 (232)  8 19  2  4 37.6   42.105    73.600   19  6.5 1992232020437
...
*/

#include"header.h"

Event *read_event_file (char *filename, int *n_events){

  FILE *fp;
  char line[LINE_LEN];
  int i, list_length;
  Event *event_list;

  fp=fopen(filename, "r");

  // Read the first line (header of file)
  fgets(line, LINE_LEN, fp);

  // Allocate event list to be N_EVENTS long
  list_length=N_EVENTS;
  event_list=malloc(list_length*sizeof(Event));

  // Read each subsequent line up to the end of the file and
  // add to event list

  i=0;
  while(fgets(line, LINE_LEN, fp)!=NULL){
 
    if(i == list_length){
      // have reached limit of event_list
      list_length*=2;
      event_list=realloc(event_list, list_length*sizeof(Event));
    }

    sscanf(line,"%d (%d) %d %d %d %d %f %f %f %f %f %s",&event_list[i].year, &event_list[i].jday, &event_list[i].month, &event_list[i].day, &event_list[i].hour, &event_list[i].min, &event_list[i].sec, &event_list[i].lat, &event_list[i].lon, &event_list[i].depth, &event_list[i].mag, event_list[i].identifier);

    i++;
  }

  *n_events=i;
  event_list=realloc(event_list, *n_events*sizeof(Event));

  fclose(fp);

  return event_list;
}
