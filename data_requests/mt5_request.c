#include"header.h"

int mt5_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user){

  int i,j, n_requests;
  char filename[REQ_NAME_LEN];
  time_t origin_time, end_time;
  
  for (i=0;i<n_events;i++){

    // Create breqfast filename
    sprintf(filename,"%s.BREQfast",event_list[i].identifier);

    // Add header information to breqfast file
    add_breqfast_header(filename, user);

    // Create start time and stop time of request
    origin_time=get_origin_time_t(event_list[i]);
    end_time=add_to_time_t(origin_time, MT5_SECONDS);

    // Search through the station list for the stations 
    n_requests=0;

    for(j=0;j<n_stations;j++){

      if(satisfy_mt5_condition(station_list[j], event_list[i])){

        write_request(filename, station_list[j].station, station_list[j].network, &origin_time, &end_time, 3, "BHZ BHE BHN");

        n_requests++;

      }

    }

    printf("File %s contains %d request lines\n", filename, n_requests);
  }
}

int satisfy_mt5_condition(Station station, Event event){

   float elat, elon, slat, slon, azm, bzm, ddg, dkm;
   int satisfaction;

   elat=event.lat;
   elon=event.lon;
   slat=station.lat;
   slon=station.lon;
   
   distaz_(&elat, &elon, &slat, &slon, &azm, &bzm, &ddg, &dkm);

   satisfaction=0;

   if(ddg <= MT5_S_MAX && ddg >= MT5_S_MIN)
     satisfaction=1;

   return satisfaction;
}

int geoscope_mt5_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user){

  int i,j, n_requests;
  char filename[REQ_NAME_LEN];
  time_t origin_time, end_time;
  
  for (i=0;i<n_events;i++){

    // Create netdc filename
    sprintf(filename,"%s.NETDC",event_list[i].identifier);

    // Add header information to netdc file
    add_breqfast_header(filename, user);

    // Create start time and stop time of request
    origin_time=get_origin_time_t(event_list[i]);
    end_time=add_to_time_t(origin_time, MT5_SECONDS);

    // Search through the station list for the stations 
    n_requests=0;

    for(j=0;j<n_stations;j++){

      if(satisfy_mt5_condition(station_list[j], event_list[i])){

        write_geoscope_request(filename, station_list[j].station, &origin_time, &end_time, "BHZ BHE BHN");

        n_requests++;

      }

    }

    printf("File %s contains %d request lines\n", filename, n_requests);
  }
}
