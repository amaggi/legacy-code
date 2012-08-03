#include"header.h"

int range_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user){

  int i,j, n_requests, n_seconds, n_channels;
  char filename[REQ_NAME_LEN], channel_line[LINE_LEN], channel[4];
  time_t origin_time, end_time;
  float min_range, max_range;

  printf("Give minimum and maximum range in degrees: ");
  scanf("%f %f", &min_range, &max_range);
  printf("Give time legth of seismogram in seconds (integer): ");
  scanf("%d", &n_seconds);
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
  
  
  for (i=0;i<n_events;i++){

    // Create breqfast filename
    sprintf(filename,"%s.BREQfast",event_list[i].identifier);

    // Add header information to breqfast file
    add_breqfast_header(filename, user);

    // Create start time and stop time of request
    origin_time=get_origin_time_t(event_list[i]);
    end_time=add_to_time_t(origin_time, n_seconds);

    // Search through the station list for the stations 
    n_requests=0;

    for(j=0;j<n_stations;j++){

      if(satisfy_range_condition(station_list[j], event_list[i], min_range, max_range)){

        write_request(filename, station_list[j].station, station_list[j].network, &origin_time, &end_time, n_channels, channel_line);

        n_requests++;

      }

    }

    printf("File %s contains %d request lines\n", filename, n_requests);
  }
}

int satisfy_range_condition(Station station, Event event, float min_range, float max_range){

   float elat, elon, slat, slon, azm, bzm, ddg, dkm;
   int satisfaction;

   elat=event.lat;
   elon=event.lon;
   slat=station.lat;
   slon=station.lon;
   
   distaz_(&elat, &elon, &slat, &slon, &azm, &bzm, &ddg, &dkm);

   satisfaction=0;

   if(ddg <= max_range && ddg >= min_range)
     satisfaction=1;

   return satisfaction;
}

int geoscope_range_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user){

  int i,j, n_requests, n_seconds, n_channels;
  char filename[REQ_NAME_LEN], channel_line[LINE_LEN], channel[4];
  time_t origin_time, end_time;
  float min_range, max_range;

  printf("Give minimum and maximum range in degrees: ");
  scanf("%f %f", &min_range, &max_range);
  printf("Give time legth of seismogram in seconds (integer): ");
  scanf("%d", &n_seconds);
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
  
  
  for (i=0;i<n_events;i++){

    // Create breqfast filename
    sprintf(filename,"%s.NETDC",event_list[i].identifier);

    // Add header information to breqfast file
    add_breqfast_header(filename, user);

    // Create start time and stop time of request
    origin_time=get_origin_time_t(event_list[i]);
    end_time=add_to_time_t(origin_time, n_seconds);

    // Search through the station list for the stations 
    n_requests=0;

    for(j=0;j<n_stations;j++){

      if(satisfy_range_condition(station_list[j], event_list[i], min_range, max_range)){

        write_geoscope_request(filename, station_list[j].station, &origin_time, &end_time, channel_line);

        n_requests++;

      }

    }

    printf("File %s contains %d request lines\n", filename, n_requests);
  }
}
