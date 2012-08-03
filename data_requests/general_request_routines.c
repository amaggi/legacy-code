#include"header.h"

time_t get_origin_time_t(Event event){

  time_t origin;
  struct tm origin_tm;

  origin_tm.tm_sec = (int)event.sec; 
  origin_tm.tm_min = event.min;
  origin_tm.tm_hour = event.hour;
  origin_tm.tm_mday = event.day;
  origin_tm.tm_mon = event.month-1;
  origin_tm.tm_year = event.year-1900;

  origin=mktime(&origin_tm);

  return origin;

}

time_t add_to_time_t(time_t starttime, int n_seconds){

  time_t endtime;

  endtime = starttime + (long) n_seconds;

  //printf("%ld + %d = %ld\n",starttime, n_seconds, endtime);

  return endtime;
}

void add_breqfast_header(char *filename, char *user){

// The breqfast profile file user contains all the header information
// other than Label and .END which are added here
/* E.g.
.NAME Alessia Maggi
.INST University of Cambridge 
.MAIL Bullard Labs., Madingley Rise, Madingley Rd, Cambridge CB3 0EZ, UK 
.EMAIL am208@cam.ac.uk
.PHONE +44-1223-337182 
.FAX +44-1223-360779
.MEDIA FTP 
.ALTERNATE MEDIA EXABYTE - 2 gigabyte
*/


  FILE *fp, *fpin;
  char line[LINE_LEN];
  
  fpin=fopen(user, "r");
  fp=fopen(filename, "w");

  while(fgets(line, LINE_LEN, fpin)!=NULL)
    fprintf(fp,"%s", line);

  fprintf(fp,".LABEL SEED.%s\n", filename);
  fprintf(fp,".END\n");

  fclose(fpin);
  fclose(fp);
 
}

void write_request(char *filename, char *station, char *network, time_t *origin, time_t *endtime, int n_cha, char *channel_string){

// Format description of request file from the IRIS website
/*
STA NN YYYY MM DD HH MM SS.T  YYYY MM DD HH MM SS.T #_CH CH1 CH2 CHn

where:

 STA  is the station
 NN   is network code 
 YYYY is the year - 2 digit entries will be rejected!
 MM   is the month
 DD   is the day
 HH   is the hour
 MM   is the minute
 SS.T is the second and tenths of seconds
 #_CH is the number of channel designators in the immediately following list
 CHn  is a channel designator that can contain wildcards
*/

  FILE *fp;
  char otime[TIME_LEN], etime[TIME_LEN];
  struct tm *otime_tm, *etime_tm;

  fp=fopen(filename, "a");

  strftime(otime, (size_t) TIME_LEN, "%Y %m %d %H %M %S.0" , localtime(origin));
  strftime(etime, (size_t) TIME_LEN, "%Y %m %d %H %M %S.0" , localtime(endtime));

  fprintf(fp,"%-4.4s %-2.2s %s %s %d %s\n", station, network, otime, etime, n_cha, channel_string);

  fclose(fp);

}

void write_geoscope_request(char *filename, char *station, time_t *origin, time_t *endtime, char *channel_string){

// Format description of request file from the geoscope
/*
  Une ou plusieurs lignes de requete(une requete par ligne):
- pour une demande d'inventaire
.INV <DATA_CENTER> <NETWORK> <STATION> <LOCATION> <CHANNELS> <START_TIME> <END_TIME>
- pour une demande de reponses instrumentales
.RESP <DATA_CENTER> <NETWORK> <STATION> <LOCATION> <CHANNELS> <START_TIME> <END_TIME>
- pour une demande de donnees
.DATA <DATA_CENTER> <NETWORK> <STATION> <LOCATION> <CHANNELS> <START_TIME> <END_TIME>

avec
<DATA_CENTER> le centre de donnees distribuant les donnees d'un ou plusieurs
reseaux. Les centres de donnees connus a ce jour: GEOFON GEOSCOPE, IRIS_DMC,
IRIS_DMC2, NCEDC (Berkeley), ORFEUS

<NETWORK> le reseau (code reseau du SEED, cf. annexe J). Les reseaux connus a
ce jour: AK, AS, AU, AZ, BK, CD, CI, CN, CS, DW, GE, G, GN, GR, HG, II, IU, KN,
MN, MX, NC, NR, PS, RS, SE, SR, TS, US, UU, UW

<STATION> la station ou liste de stations (separateur d'element de liste:
blanc, delimiteur de liste: '  ').  Exemple de stations a Geoscope (cf leWEB
geoscope.ipgp.jussieu.fr) : KIP, SSB, ATD, ...

<LOCATION> numero du numeriseur de la station. Exemple: *, 01, 02

<CHANNEL> le canal ou la liste de canal. Le nommage des canaux respecte la convention du SEED (cf.  annexe A).  Exemple: BHZ, BHN, BHE ou si les composantes nord et ne sont pas correctement orientees BHZ, BH2, BH3

<START_TIME> et <END_TIME> heures de debut et de fin delimitant la periode
voulue. Grammaire d'une date: "YYYY MM DD HH MM SS.FFFF" (ne pas oublier les
quotes!).  Exemple: "1995 03 24 22 40 00" pour 22:40 min et 00 seconde le
24/03/1995 

----------------------------------------------------

Recommandation:
1> respectez les quotes dans les lignes de requete(elles definissent les arguments de la requete comme une chaine de caractere)
*/

  FILE *fp;
  char otime[TIME_LEN], etime[TIME_LEN];
  struct tm *otime_tm, *etime_tm;

  fp=fopen(filename, "a");

  strftime(otime, (size_t) TIME_LEN, "%Y %m %d %H %M %S" , localtime(origin));
  strftime(etime, (size_t) TIME_LEN, "%Y %m %d %H %M %S" , localtime(endtime));

  fprintf(fp,".DATA GEOSCOPE G %-4.4s * \"%s\" \"%s\" \"%s\"\n", station, channel_string, otime, etime);

  fclose(fp);

}
