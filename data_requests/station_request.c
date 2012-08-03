#include"header.h"

int station_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user){

  int i,j, n_requests, timelength, n_channels;
  char filename[REQ_NAME_LEN], network[NET_LEN];
  char station[STATION_LEN], channel[4], channel_line[LINE_LEN];
  time_t origin_time, end_time;

  printf("Give network name:\n");
  scanf("%s", network);
  printf("Give station name:\n");
  scanf("%s", station);
  printf("Give time in seconds:\n");
  scanf("%d", &timelength);
  printf("Give number of channels required: ");
  scanf("%d", &n_channels);
  strcpy(channel_line,"");
  for (i=0;i<n_channels;i++){
    printf("Give channel %d : ", i+1);
    scanf("%s", channel);
    strcat(channel_line, channel);
    strcat(channel_line, " ");
    fprintf(stderr,"%s\n", channel_line);
  }


  // Create breqfast filename
  sprintf(filename,"%s-%s.BREQfast",network,station);

  // Add header information to breqfast file
  add_breqfast_header(filename, user);


    n_requests=0;
  
  for (i=0;i<n_events;i++){


    // Create start time and stop time of request
    origin_time=get_origin_time_t(event_list[i]);
    end_time=add_to_time_t(origin_time, timelength);

    // Search through the station list for the stations 

    for(j=0;j<n_stations;j++){

      if(strncmp(network, station_list[j].network, NET_LEN)==0 && strncmp(station, station_list[j].station, STATION_LEN)==0){

        write_request(filename, station_list[j].station, station_list[j].network, &origin_time, &end_time, n_channels, channel_line);

        n_requests++;

      }

    }

  }
    printf("File %s contains %d request lines\n", filename, n_requests);
}

int geoscope_station_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user){

  int i,j, n_requests, timelength, n_channels;
  char filename[REQ_NAME_LEN], network[NET_LEN];
  char station[STATION_LEN], channel[4], channel_line[LINE_LEN];
  time_t origin_time, end_time;

  printf("Give network name:\n");
  scanf("%s", network);
  printf("Give station name:\n");
  scanf("%s", station);
  printf("Give time in seconds:\n");
  scanf("%d", &timelength);
  printf("Give number of channels required: ");
  scanf("%d", &n_channels);
  strcpy(channel_line,"");
  for (i=0;i<n_channels;i++){
    printf("Give channel %d : ", i+1);
    scanf("%s", channel);
    strcat(channel_line, channel);
    strcat(channel_line, " ");
    fprintf(stderr,"%s\n", channel_line);
  }


  // Create breqfast filename
  sprintf(filename,"%s-%s.NETDC",network,station);

  // Add header information to breqfast file
  add_breqfast_header(filename, user);


    n_requests=0;
  
  for (i=0;i<n_events;i++){


    // Create start time and stop time of request
    origin_time=get_origin_time_t(event_list[i]);
    end_time=add_to_time_t(origin_time, timelength);

    // Search through the station list for the stations 

    for(j=0;j<n_stations;j++){

      if(strncmp(network, station_list[j].network, NET_LEN)==0 && strncmp(station, station_list[j].station, STATION_LEN)==0){

        write_geoscope_request(filename, station_list[j].station, &origin_time, &end_time, channel_line);

        n_requests++;

      }

    }

  }
    printf("File %s contains %d request lines\n", filename, n_requests);
}
