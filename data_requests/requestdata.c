// Program to request data from Iris - essentially to generate breqfast files

#include"header.h"

int main(int argc, char *argv[]){

  int n_stations, n_events, selection_type;
  Station *station_list;
  Event *event_list;

  if(argc != 4){
    printf("Usage requestdata CMTSOLUTION fsdn.dat userprofile\n");
    return;
  }

//  event_list=read_event_file(argv[1],&n_events);
//  station_list=read_station_time_file(argv[2],&n_stations);
  event_list=read_CMTSOLUTION_file(argv[1],&n_events);
//  station_list=read_STATIONS_file(argv[2],&n_stations);
  station_list=read_station_file(argv[2],&n_stations);

  printf("Have read %d stations and %d events from files.\n", n_stations, n_events);

  printf("Select type of request:\n");
  printf("%d\tMt5 selection\n",MT5);
  printf("%d\tStation selection\n",STATION);
  printf("%d\tDistance range selection\n",RANGE);
  scanf("%d",&selection_type);

  switch(selection_type){
    case(MT5):
      mt5_request(station_list, n_stations, event_list, n_events, argv[3]);
      break;
    case(STATION):
      station_request(station_list, n_stations, event_list, n_events, argv[3]);
      break;
    case(RANGE):
      range_request(station_list, n_stations, event_list, n_events, argv[3]);
      break;
    default:
      printf("Incorrect selection type\n");
      return;
      break;
  }


  free(event_list);
  free(station_list);

}
