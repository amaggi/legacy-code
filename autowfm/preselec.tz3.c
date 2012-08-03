#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
void readper();
int tt[6];
main()
{
/* 7 juin 2000 : version simplifiee de preselec pour le
   traitement des Rayleigh seuls */
int ndata,i,iray,ilovray,ista,itest;
char * bido;
char * bidofin;
char ligne[100];
char com[100];
char com2[200];
char com3[200];
char mouvons[100];
char data[100];
char dataname[100];
char nomin1[100];
char nomin2[100];
char nomso1[100];
char nomso4[100];
FILE *in1;
FILE *so1;
FILE *so4;

/* Ouverture des 2 fichiers de sortie
   contenant les valeurs de filtrage des 
   intercorelogrammes */
strcpy(nomso1,"DATA.TFILRAYray");
if(!(so1=fopen(nomso1,"w")))
       {printf ("\n pb opening file : %s",nomso1);
        exit(1);}

strcpy(nomso4,"DATA.TFILRATES");
if(!(so4=fopen(nomso4,"w")))
       {printf ("\n pb opening file : %s",nomso4);
        exit(1);}

/* Ouverture du fichier d'entree */
strcpy(nomin1,"DATARAYray.LST");
if(!(in1=fopen(nomin1,"rt")))
       {printf ("\n pb opening file : %s",nomin1);
        exit(1);}

/*fgets(ligne,100,in1);
sscanf(ligne,"%d ",&ndata);*/
 
 fscanf(in1,"%d",&ndata);
       
for(i=0;i<ndata;i++)
  {
  iray=0;
/*fgets(ligne,100,in1); */
  fscanf(in1,"%s",ligne);
  printf("\n %s\n",ligne);

  strcpy(dataname,ligne);

/*bidofin=ligne + strcspn(ligne,"SAC")+3;
  strcpy(bidofin,"");
  strcpy(dataname,bido);
  printf("\n dataname : %s strcspn(ligne,SAC)+3 %d",dataname,strcspn(ligne,"SAC")+3);*/

  strcpy(com2,"newsac ray1.50.m dat ");
  strcpy(com3,"newsac ray2.50.m dat ");
  strcpy(com,dataname);
  strcat(com2,com);
  strcat(com3,com);
  printf("\n %s \n",com2);
  system(com2);                                  /*lancement macro  ray */
  system("rm bid1.sac bid2.sac bid3.sac");  
  printf("\n %s \n",com3);
  system(com3);                                  /*lancement macro  ray */
  system("rm bid1.sac bid2.sac bid3.sac");  
  strcpy(nomin2,"vfilt");
  readper(nomin2);
  system("rm vfilt");
/* NB : tt[0]=40s (harmoniques)
        tt[1]=50s;tt[2]=70s;tt[3]=90s     */
  if(((tt[1]==1) && (tt[2]==1)) || ((tt[1]==1) && (tt[3]==1)))
      {
      iray=1;
      fprintf (so1,"\n  %s %d %d %d %d %d %d",dataname,tt[0],tt[1],tt[2],tt[3],tt[4],tt[5]);
/*    don't want to copy data  (sylvana)
      strcpy(mouvons,"cp ");
      strcat(mouvons,dataname);
      strcat(mouvons," /export/home/pilidou/TMP/test");
      printf("\n  MOUVONS iray=1 :  %s \n",mouvons);
      system(mouvons);   */
      }
      else
      {
      iray=0;
      fprintf (so4,"\n %s %d %d %d %d %d %d",dataname,tt[0],tt[1],tt[2],tt[3],tt[4],tt[5]);
      }
  }
fclose (in1);
fclose (so1);
}
/*#####reading the frequency file############*/
void readper(nomin2)
char nomin2[100];
{
float fr;
char ligne[100];
char bid[100];
char bid2[100];
FILE *in2;
int ij;
if(!(in2=fopen(nomin2,"rt")))   
       {printf ("\n pb opening file vf %s \n",nomin2);
        exit(1);}
  for(ij=0;ij<6;ij++)   
       {  
      /*fgets(ligne,100,in2);
        sscanf(ligne,"%s %s %e",bid,bid2,&fr);*/
        fscanf(in2,"%*s %*s %e",&fr);
        tt[ij]=(int)fr;
       }  
fclose (in2);
}
