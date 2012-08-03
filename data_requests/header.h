#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>

#define MT5 1
#define STATION 2
#define RANGE 3

#define N_STATIONS 100
#define N_EVENTS 100

#define NET_LEN 2
#define STATION_LEN 5
#define LOCATION_LEN 120
#define LINE_LEN 255
#define TIME_LEN 22

#define EV_IDENT_LEN 13
#define REQ_NAME_LEN 22 // EV_IDENT_LEN+".BREQfast"

#define MT5_P_MIN 30
#define MT5_S_MIN 35
#define MT5_P_MAX 90
#define MT5_S_MAX 85

#define MT5_SECONDS 1800

typedef struct STation {
  char network[NET_LEN];
  char station[STATION_LEN];
  float lat, lon, elevation, depth;
  char location[LOCATION_LEN];
} Station;

typedef struct EVent {
  int year, jday, month, day, hour, min;
  float sec, lat, lon, depth, mag;
  char identifier[EV_IDENT_LEN];
} Event;

Station * read_station_file (char *filename, int *n_stations);
Station * read_station_time_file (char *filename, int *n_stations);
Event * read_event_file (char *filename, int *n_events);
Station * read_STATIONS_file (char *filename, int *n_stations);
Event * read_CMTSOLUTION_file (char *filename, int *n_events);

int mt5_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user);
int geoscope_mt5_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user);
int station_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user);
int geoscope_station_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user);
int range_request(Station station_list[], int n_stations, Event event_list[], int n_events, char *user);

int satisfy_mt5_condition(Station station, Event event);
int satisfy_range_condition(Station station, Event event, float min_range, float max_range);

void distaz_(float *elat, float *elon, float *slat, float *slon, float *azm, float *bzm, float *ddg, float *dkm);

time_t get_origin_time_t(Event event);
time_t add_to_time_t(time_t starttime, int n_seconds);

void add_breqfast_header(char *filename, char *user);
void write_request(char *filename, char *station, char *network, time_t *origin, time_t *endtime, int n_cha, char *channel_string);
void write_geoscope_request(char *filename, char *station, time_t *origin, time_t *endtime, char *channel_string);

