#include <stdio.h> /*include standard input/output library */
#include <math.h> /* include math library */
#include <malloc.h>
#include <string.h>  /*include string handling utilities*/
#include <stdlib.h>

/* file from which to read paths from */
#define pathfile "WFMPATHS"


/*######################## declaration of structure parameters #################################*/
struct cmt    
{
int axepp;
int axepa;
int axetp;
int axeta;
float evdp;
double moment;
float dist;
float azim;
float dur;
};
/*####################### declaration of functions ##################################*/
void initialisation();
void initperiode();
void readcmt();
void initnewrun();
void classement();
float azimuth();
int residus();
/*############################ global variables ###################################*/
#define Epsilon 0.00000001
#define pi 3.14159265358979323846264338
int tpmaxray;
int tpmaxlov;
int tpminray;
int tpminlov;
int tmax;
int tmin;
/*############################ Comments #####################################*/
/*
  Code by Eric Debayle, changed by Sylvana Pilidou, March 2002.
  

  VERSION Fevrier 98 retenue pour la modelisation des Rayleigh seuls
- les periodes>100s ne sont pas inversees pour SKIPPY.
- La contribution de chaque mode par rapport a l'energie
  totale du synthetique est utilisee pour selectionner
  les enveloppes.
  
== he contribution of each mode compared to the total energy of synthetic is utilisee for selectionner the envelopes
  
  - Les enveloppes sont inversees a longue, puis cp, puis on
  passe a la phase lp et cp. 2 lobes peuvent etre selectiones
  sur chaque enveloppe l'erreur sur la donnees etant dependante
  du rapport signal/bruit.

  == the envelopes are inversees has long, then CP, then one passes has the phase LP and CP 2 lobes can be selectiones on each envelope the error on the data being dependante signal-to-noise report/ratio.

- Le dephasage entre intercore synthe et reels est teste.
  La valeur minimale du dephasage est signalee pour le mode
  et la periode lorsque celui-ci est >pi. Pourra etre
  utilise + tard pour trier des donnees.

== The dephasage between intercorell of the  synth. and real is tested. The minimal value of the dephasage is sign for the mode and the period when this one is > pi. Could be uses + late for sorting data.

Choix des criteres d'acceptation
- Stabilite du modele calculee  partir de khim
- enveloppes : 90% redvar ou qui<3 que ce soit lp ou cp
- reduction d'energie: soit>90% soit moyenne Love/Rayleigh>90%

== Choice of the criteria of acceptance 
   - Stability of the model calculee to start from khim 
    - envelopes: 90% redvar or qui<3 that it is LP or CP 
    - reduction of energy: soit>90% are average Love/Rayleigh>90% 

NB: dans cette version, readcmt a ete modifie pour lire
    les fichiers CMT harvard dont le format a legerement change
    apres juillet 96.

== NB: This version reads new CMT data format (changed in 1996)
	
############################ MAIN PROGRAM  ########################################*/
void main()
{
char ligne[100];
char adrdatray[80],adrdatlov[80];
char nomdatray[45],nomdatlov[45],savetv[45];
char chin[12],chite[9],chout[8],chorder[30],bids[20];
char yiray[200],yilov[150],souray[200],soulov[150],dpray[200],dplov[200];
int ilr,ipath,npaths,idata,irun;
int itstart,itnewstart,iteend,icompt,npts;
int periodr[7],periodl[7];
float pas,dista;
FILE * dataray;
FILE * datalov;
FILE * fichdat;

char path_data[150], path_cmt[150], path_saito[150], path_resp[150] ;
FILE *fd ;

/* sylvana: open paths file and get paths */
  if (!(fd = fopen(pathfile,"r"))) {
     printf ("Cannot open file WFMPATHS\n");
     exit(1);
  }
fscanf(fd, "%s", path_data);/*read first line of file fd*/
fscanf(fd, "%s", path_cmt);/*read second line of file fd*/
fscanf(fd, "%s", path_saito);/*read third line of file fd*/
fscanf(fd, "%s", path_resp);/*read fourth line of file fd*/

fclose(fd);

printf("path_data contains %s\n", path_data);
printf("path_cmt contains %s\n", path_cmt);
printf("path_saito contains %s\n", path_saito);
printf("path_resp contains %s\n", path_resp);

/* Nouvelle version : inversion de la phase a LP. 2 Phases
   sont introduites dans l'inversion.
   idata est un indicateur sur le type de donnees inversees
   qd idata=-2  on inverse seulement les enveloppes LP
   qd idata=-1  on inverse les enveloppes LP + CP         
   qd idata=0  on inverse les enveloppes + phase LP            
   qd idata=+1 on inverse les enveloppes + phase CP           

== This version: inversion of phase of LP. 2 phases instroduced in inversion.
   idata indicated the type of inversion
   idata = -2  invert only envelopes LP
   idata = -1  invert envelopes LP + phase CP
   idata = 0   invert envelopes + phase LP
   idata = 1   invert envelopes + phase CP
*/


printf ("\n  If we work with more than 5 modes, nmode has to be");
printf ("\n  initialised in initperiode");
printf ("\n  We work with Love and Rayleigh (1) or only Rayleigh(0) \t");
gets(ligne);   
sscanf(ligne,"%d ",&ilr);
/* if Rayleigh, check that  DATATOTRAY exists/readable, otherwise stop*/
if(ilr==0) 
      {  
        if(!(dataray=fopen("DATATOTRAY","r")))
           {printf ("\n cannot open file DATATOTRAY");
            exit(1);}
        strcpy(nomdatlov,"bid");
      }  
/* if Rayleigh, check that  DATATOTRAY exists/readable, otherwise stop*/
else if (ilr==1) 
      {  
        if(!(dataray=fopen("DATATOTRAY","r")))
           {printf ("\n cannot open file DATATOTRAY");
            exit(1);}  
        if(!(datalov=fopen("DATATOTLOV","r")))   
           {printf ("\n cannot open file DATATOTLOV"); 
            exit(1);}      
      } 

/* Number of seismograms we have ready (with dp and yi files */
printf ("\n number of data for which we have :");
printf ("\n - a .dat file (real seismogram)");
printf ("\n - a dp and yir files for both the average");
printf ("\n   structure and the source region\t");
gets(ligne);
sscanf(ligne,"%d ",&npaths);
/* ****************************************************/
/* loop on seismograms and read the DATATOTRAY file */
/* ****************************************************/
for (ipath=1;ipath<=npaths;ipath++)
   {
   tpmaxray=0;
   tpminray=0;
   tpmaxlov=0;
   tpminlov=0;
   idata=-2;
   fgets(ligne,100,dataray);
   sscanf(ligne,"%s %d %d %d %d %d %d",nomdatray,&periodr[0],&periodr[1],
          &periodr[2],&periodr[3],&periodr[4],&periodr[5]);
   
   sprintf(adrdatray, "%s%s",path_data,nomdatray);
   if(!(fichdat=fopen(adrdatray,"r")))
       {printf ("\n cannot open file: %s%s",path_data,nomdatray);
        exit(1);}
/* One will seek the distance which avere necessary for the choice of the periodes in initperiode. Fichdat is read again initialization. One could simplify the programming here by reading fichdat only one time. To avoid the pbs I prefere leaving like Ca for the moment*/
   fgets(ligne,100,fichdat);  /* On va chercher la distance qui s'avere              */
   fgets(ligne,100,fichdat);  /* necessaire pour le choix des periodes dans          */
   fgets(ligne,100,fichdat);  /* initperiode. Fichdat est relu dans initialisation   */
   fgets(ligne,100,fichdat);  /* On pourrait simplifier le programmation en ne lisant*/
   fgets(ligne,100,fichdat);  /* fichdat qu'une fois ici. Pour eviter les pbs je     */
   fgets(ligne,100,fichdat);  /* prefere laisser comme ca pour le moment.            */
   fgets(ligne,100,fichdat);
   sscanf(ligne,"%d %f %f %s",&npts,&pas,&dista,bids);
   fclose (fichdat);

/* SUBROUTINE CALL ************************************************  */
   initperiode(dista,periodr,idata,0);
/******************************************************************* */   
   tpmaxray=tmax;
   tpminray=tmin;
/* also read love waves file if have both */   
   if(ilr==1)
        {fgets(ligne,100,datalov);
         sscanf(ligne,"%s %d %d %d %d %d %d",nomdatlov,&periodl[0],&periodl[1],
         &periodl[2],&periodl[3],&periodl[4],&periodl[5]); 
         initperiode(dista,periodl,idata,1);
         tpmaxlov=tmax;
         tpminlov=tmin;}
   sprintf(adrdatlov,"%s%s",path_data,nomdatlov);

   itstart=0;
   itnewstart=0;
   iteend=3;
   icompt=1; 
   initialisation(path_data,path_cmt,path_saito,path_resp,nomdatray,nomdatlov,adrdatray,adrdatlov,ilr,yiray,yilov,souray,soulov,dpray,dplov);
          
   printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
   printf("\n             TRAJET no %d on %d ",ipath,npaths);
   printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
//   system ("rm autorun");
   system ("cat initrun corprun>autorun");
   system ("chmod 744 autorun");

   strcpy(chin,"./autorun");
   debut:
   strcpy(chorder,chin);
   sprintf(chite," %d %d > ",itstart,iteend);
   strcat(chorder,chite);
   sprintf(chout,"error%d",icompt);
   strcat(chorder,chout);

   printf("\n---------------------------------------");
   printf("\n %s ",chorder);
   printf("\n---------------------------------------\n");

   system (chorder); 

   irun=residus(chout,nomdatray,itnewstart+1,iteend+1,idata,ilr);

   printf ("\n ---------------------------------");
   printf ("\n FIN run %d IRUN = %d",icompt,irun);
   printf ("\n ---------------------------------\n");

/* on sort: rapport synthe/reel pas bon 
 == synthe/reel ratio not good  */
   if (irun<0)                         
     { goto sortie;}

/* 
on continue : pas assez d'ite pour expliquer les donnees
ou modele pas stabilise 
== not enough of ... to explain the data or model not stable  */
   if (irun==0) 
   { 
     if( ((iteend<13) && (idata<=0.)) ||   
         ((iteend<16) && (idata>0.)) )
       { 
        itstart=iteend+1;                                               
        iteend=itstart+2; 
        icompt++;
        if(idata>=0.) strcpy(chin,"./autonewrun");
        if(idata<0.)strcpy(chin,"./autorun");
        goto debut;
       }
     else 
        goto sortie;     
    }              
                                         
/*
irun>0
cas 1: on a inverse les env lp seulement      => on rajoute les enveloppes cp 
cas 2: on a inverse les env lp +cp            => on rajoute la phase lp      
cas 3: on a inverse toutes les env+phase lp   => on rajoute la phase cp
cas 4: on a deja inverse enveloppes +phase    => on sort             
*/
   if (irun>0)                             
   {  
     if (idata<=0.)
      {  itstart=iteend+1;
         itnewstart=itstart;
         iteend=itstart+2;  
         icompt++;        
         initnewrun(itstart);
         system("cat initrun newcorp>autonewrun");
         system("chmod 744 autonewrun");
         if (idata==-2)
             { system("mv pick.ray.1 pick.ray.1.lp");
               if(ilr==1)system("mv pick.lov.1 pick.lov.1.lp");}
         if (idata==-1)
             { system("mv pick.ray.1 pick.ray.1.env");
               if(ilr==1)system("mv pick.lov.1 pick.lov.1.env");}
         if (idata==0)
             { system("mv pick.ray.1 pick.ray.1.plp");
               if(ilr==1)system("mv pick.lov.1 pick.lov.1.plp");}
         idata++;
         strcpy(chin,"./autonewrun");
         sprintf(savetv,"cp tv.it%d tv.it%ds",itstart,itstart);
         system(savetv);          
         printf("\n SAVETV  %s",savetv);
         initperiode(dista,periodr,idata,0);
         if(ilr==1) initperiode(dista,periodl,idata,1);
         goto debut;                          
       }
      else if (idata>0)                
         goto sortie;         
       }
   system ("rm bid");
   sortie:;
   classement(irun,nomdatray,iteend,yiray,yilov,souray,soulov,dpray,dplov);
 }
fclose(dataray);
if (ilr==1) fclose(datalov); 
}

/****************************************************/
/* end loop on seismograms and read the DATATOTRAY file */
/****************************************************/
/*############################ END MAIN PROGRAM  ########################################*/

/*############################ SUBROUTINES #############################################*/
int residus(chin,nomdatray,itstart,iteend,idata,ilovray)
int iteend, itstart,ilovray,idata;
char * nomdatray;
char * chin;
{
float quipha;
float energie,energiein,redenergie,avenergie,averedenergie;
float echs[4], qui[20],khim[20],quiinit[20], nergie[2][20];
float quienv[20],quienvi[20];
double redvar,redvarlp;
int i,j,irun,il,ienergie;
char chout[45];
char bids[120];
char ligne[120];
char commande[300];
FILE * out;
FILE * in;

/* Routine d'evaluation du resultat de la wfm inversion.
   Ne marche que si on inverse un seul Rayleigh ou en anisotrope
   un seul love et un seul rayleigh soit 2 sismogrammes.  */

if(!(in=fopen(chin,"r")))
    {printf ("\n routine residus pb opening file %s :",chin);
     exit(1);}
il=strlen(nomdatray);
if (il<=4) 
    {printf ("\n Le nom du fichier dat est mal interprete!!!!");
     exit(1);}
strcpy(chout,"res.");
strncat(chout,nomdatray,il-4);
sprintf(commande,"rm %s",chout);
system (commande);
if(!(out=fopen(chout,"w")))
    {printf ("\n routine residus pb opening file %s :",chout);
     exit(1);}
redvar=0.;
redvarlp=0.;
redenergie=0.;
irun=0;

/*lecture facteur d'echelle*/
j=1;
do
   fgets(ligne,120,in);
while(strstr(ligne,"facteur d'echelle")==0);
fputs(ligne,out);
if (ilovray==0)
     sscanf(ligne,"%s %s %f",bids,bids,&echs[j]); 
if (ilovray==1)
     sscanf(ligne,"%s %s %s %f",bids,bids,bids,&echs[j]);
printf("\n FACT %f",echs[j]);
if (ilovray==1)
{j=2;
 fgets(ligne,120,in);
 fputs(ligne,out);
 sscanf(ligne,"%s %s %s %f",bids,bids,bids,&echs[j]);
 printf("\n FACT %f",echs[j]);}

for(j=1;j<=ilovray+1;j++)
    {if (echs[j]==9999.)
         {irun=-1;
         fgets(ligne,120,in);
         fputs(ligne,out);
         fprintf(out,"\nIRUN=-1 :arret imprevu du pgm\n");
         goto sortie;     } 

     if ((echs[j]>5.) || (echs[j]<0.2)) 
         {irun=-1;
          fprintf(out,"\nIRUN=-1 Valeur de FACT excessive\n");
          goto sortie;}    }

/*lect degres de libertes */
for (i=1;i<=iteend;i++)         
   { fgets(ligne,120,in);
     fputs(ligne,out); }

/*lect du qui global et khim*/
for (i=1;i<=iteend;i++)
   { fgets(ligne,120,in);      
      sscanf(ligne,"%s %s %s %s %s %f %s %s %f %s %f"
      ,bids,bids,bids,bids,bids,&quiinit[i],bids,bids,&qui[i],bids,&khim[i]);
    /* For tomo cf ci-dessous 
       sscanf(ligne,"%s %s %s %s %f %s %s %f %s %f"
       ,bids,bids,bids,bids,&quiinit[i],bids,bids,&qui[i],bids,&khim[i]);*/
     fputs(ligne,out);}

/*lect du qui enveloppe */
for (i=1;i<=iteend;i++)   
   { fgets(ligne,120,in);
     sscanf(ligne,"%s %s %s %s %s %s %s %f %s %s %f"
     ,bids,bids,bids,bids,bids,bids,bids,&quienvi[i],bids,bids,&quienv[i]);
  /*  For tomo cf ci-dessous 
      sscanf(ligne,"%s %s %s %s %s %s %f %s %s %f"
      ,bids,bids,bids,bids,bids,bids,&quienvi[i],bids,bids,&quienv[i]); */
     fputs(ligne,out); }

/* lect. du qui phase    */
for (i=1;i<=iteend;i++)
   { fgets(ligne,120,in);       
     sscanf(ligne,"%s %s %s %s %f",bids,bids,bids,bids,&quipha);
  /* For tomo cf ci-dessous 
     sscanf(ligne,"%s %s %s %f",bids,bids,bids,&quipha); */
     fputs(ligne,out); }
printf ("\n");
printf("\n quiinit %f qui %f quim %f quienvi %f quienv %f quipha %f"
,quiinit[iteend],qui[iteend],khim[iteend],
quienvi[iteend],quienv[iteend],quipha);

/*
la red. de variance est calculee pour chaque run du wfm
effectue a nbre de donnees constant et seulement pour les donnees
ajoutees dans l'inversion. Qd on en arrive a l'inversion
de la phase calculer une redvar n'a plus aucune sens car :
- On a deja fitte les enveloppes lp puis cp avec qui<3 ou redvar>80%
- Le meme fit d'une enveloppe peut etre obtenu avec une phase de
  l'intercorrelogramme qui change de +- pi=> si on a de la chance
  on fitte deja la phase apres l'inversion des enveloppes seule
  et en ajoutant la phase dans l'inversion, on ne changera rien au
  modele, ie on aura une reduction de variance nulle. Si on a pas de chance
  On est a +-pi avec les enveloppes seules et l'ajout de la phase
  va entrainer une grosse red. de variance.
On ne calcule pas de redvar apres inversion de la phase lp. On passera
directement a l'iteration suivante.
*/

if(idata==-2)
  {
   if(quiinit[itstart]!=0.)
   redvarlp=(pow((double)quiinit[itstart],2.) 
         -pow((double)qui[iteend],2.))/pow((double)quiinit[itstart],2.);
  }
else if (idata==-1)
  {
   if(quiinit[1]!=0.)
   redvarlp=(pow((double)quiinit[1],2.) 
         -pow((double)qui[itstart-1],2.))/pow((double)quiinit[1],2.);

   if(quienvi[itstart]!=0.)
   redvar=(pow((double)quienvi[itstart],2.)
         -pow((double)quienv[iteend],2.))/pow((double)quienvi[itstart],2.);
  }
/*si on inverse la phase lp il faut des enveloppes lp bien fittees*/
if((idata==-2) && ((qui[iteend]<3.) || (redvarlp>=0.8))) irun=1;
if((idata==-1) && ((qui[iteend]<3.) || (redvar>=0.8))) irun=1;

/* On ne prend en compte que le khim pour verifier 
la stabilite du modele
  if((idata==1) && (fabs(qui[iteend]-qui[iteend-1])<=0.05)
     && (fabs(khim[iteend]-khim[iteend-1])<=0.05)) irun=1; */

if((idata==1) && (fabs(khim[iteend]-khim[iteend-1])<=0.05)) irun=1;

/* lect. de Mo */
for (i=1;i<=iteend;i++)
   {for(j=1;j<=ilovray+1;j++)   
        { fgets(ligne,120,in);
          fputs(ligne,out); } }

/* lect de la distance de correlation */
fgets(ligne,120,in);           

/*lect de l'energie relative */
ienergie=0;
avenergie=0.;         
averedenergie=0.;         
for(j=1;j<=ilovray+1;j++)       
   {
    energiein=0.;
    energie=0.;         
    fgets(ligne,120,in);
    for (i=0;i<=iteend;i++)
      { fgets(ligne,120,in);
        sscanf(ligne,"%s %s %s %s %s %s %s %s %f",
               bids,bids,bids,bids,bids,bids,bids,bids,&nergie[j][i]);
        printf("\n ENERGIE %f",nergie[j][i]);
        fputs(ligne,out);     }
        energie=nergie[j][iteend];
        energiein=nergie[j][0];    
        printf("\n SISMO %d ENERGIE %f ENERGIEIN %f\n",j,energie,energiein);
        redenergie=(energiein-energie)/energiein;
        printf("\n REDENERGIE %f",redenergie);
        fprintf(out,"\n SISMO %d : Reduction d'energie   %f\n",j,redenergie);
        if((idata==1) && ((energie<=0.3) || (redenergie>=0.9)))  ienergie=ienergie+1;
        avenergie=avenergie+energie;
        averedenergie=averedenergie+redenergie;
                                                                            }

/*BEUG DETECTE le 19 Mars 98: la condition irun==1 n'etait
  pas dans le if suivant=>un modele pas stable pouvait etre
  accepte s'il permettait de fitter les donnees et la forme
  d'onde*/
        avenergie=avenergie/(ilovray+1);
        averedenergie=averedenergie/(ilovray+1);
        if ((ienergie==(ilovray+1)) && (idata==1) && (irun==1)) 
             irun=1;
        else if(((averedenergie>=0.9) || (avenergie<=0.3)) && (idata==1) && (irun==1))
             irun=1;
        else if(idata==1)
             irun=0;
     
/*En 2 runs on doit fitter la phase des enveloppes lp.*/
if (idata==-2) 
     fprintf(out,"\n Reduction de variance enveloppes lp %f",redvarlp);
if (idata==-1) 
     fprintf(out,"\n Reduction de variance enveloppes %f",redvar);
if (idata==0) irun=1;
fprintf(out,"\n---------------------------------------");

if((idata==-2) && ((qui[iteend]>=3.) && (redvarlp<0.8))) 
           fprintf(out,"\n IRUN=0 car qui>3 et redvarlp<0.8 apres inversion des enveloppes lp \n");
if((idata==-1) && ((qui[iteend]>3.) && (redvar<0.8))) 
           fprintf(out,"\n IRUN=0 car qui>3 et redvar<0.8 apres inversion des enveloppes \n");
/* if((idata==1) && ((fabs(qui[iteend]-qui[iteend-1])>=0.05)
                    || (fabs(khim[iteend]-khim[iteend-1])>=0.05)))  */
if((idata==1) && (fabs(khim[iteend]-khim[iteend-1])>=0.05))
           fprintf(out,"\n IRUN=0 car modele pas stabilise\n");
if (((ienergie!=(ilovray+1)) && ((averedenergie<0.9) && (avenergie>0.3))) && (idata==1))
           fprintf(out,"\nIRUN=0 car avenergie>0.3 et redenergie insuffisante\n");
if((idata<0) && (redvarlp==0.) && (redenergie==0.))
          {irun=1;
           fprintf(out,"\nIRUN=1 car pas de donnees lp selectionnees");}
fprintf(out,"\nFin residus IRUN= %d\n",irun);
fprintf(out,"\n-----------------------------------------------------------------------------\n");
sortie:;
fclose(in);
fclose(out);
return(irun);
}  
/*###################################################################################*/
void classement(irun,nomdatray,iteend,yiray,yilov,souray,soulov,dpray,dplov)
char * nomdatray;
char * yiray,* yilov,* souray,* soulov,* dpray,* dplov;
int irun,iteend;
{
int il;
char chout[41];
char commande[300];

il=strlen(nomdatray);
if (il<=4)
    {printf ("\n Le nom du fichier dat est mal interprete!!!!");
     printf ("\n dans la routine classement!!!!");
     exit(1);}
strcpy(chout,"");
strncat(chout,nomdatray,il-4);
printf("\n  CHOUT %s\n ",chout);
if(irun == 1)
 {
  sprintf(commande,"mkdir OK/%s",chout);
  system(commande);  
  sprintf(commande,"cp %s OK/%s",nomdatray,chout);
  system(commande);
  sprintf(commande,"cp des.in OK/%s",chout);
  system(commande);
  sprintf(commande,"mv autorun OK/%s",chout);
  system(commande);
  sprintf(commande,"mv autonewrun OK/%s",chout);
  system(commande);
  sprintf(commande,"mv tv.init OK/%s",chout);
  system(commande);
  sprintf(commande,"mv dpa.init* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv disp.init* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv att.init* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv sour.init* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv tf.init* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv gk.init* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv pick.* OK/%s",chout);
  system(commande);
  sprintf(commande,"mv tf.???.1.it%d OK/%s",iteend+1,chout);
  system(commande);
  sprintf(commande,"mv tf.???.1.it%d OK/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv tv.it%d OK/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv tv.it%d OK/%s",iteend-1,chout);
  system(commande);
  sprintf(commande,"mv dgk.it%d OK/%s",iteend,chout);
  system(commande); 
  sprintf(commande,"mv dgk.it%d OK/%s",iteend-1,chout);
  system(commande); 
  sprintf(commande,"mv des.it%d OK/%s/des.%s",iteend+1,chout,chout); 
  system(commande);
  sprintf(commande,"cp des.it%d OK/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv res.%s OK/%s",chout,chout);  
  system(commande);             
  sprintf(commande,"mv DEPHA* OK/%s",chout);  
  system(commande);             
  sprintf(commande,"mv energiemode* OK/%s",chout);  
  system(commande);  
  sprintf(commande,"mv error? ??.????.??.lhz OK/%s",chout);  
  system(commande);  
}
  if ((irun==0) || (irun < 0))
{
  sprintf(commande,"mkdir REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"cp %s REJECTED/%s",nomdatray,chout);
  system(commande);
  sprintf(commande,"cp des.in REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv autorun REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv autonewrun REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv tv.init REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv dpa.init* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv disp.init* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv att.init* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv sour.init* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv tf.init* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv gk.init* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv pick.* REJECTED/%s",chout);
  system(commande);
  sprintf(commande,"mv tf.???.1.it%d REJECTED/%s",iteend+1,chout);
  system(commande);
  sprintf(commande,"mv tf.???.1.it%d REJECTED/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv tv.it%d REJECTED/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv tv.it%d REJECTED/%s",iteend-1,chout);
  system(commande);
  sprintf(commande,"mv dgk.it%d REJECTED/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv dgk.it%d REJECTED/%s",iteend-1,chout);
  system(commande);
  sprintf(commande,"mv des.it%d REJECTED/%s/des.%s",iteend+1,chout,chout);
  system(commande);
  sprintf(commande,"cp des.it%d REJECTED/%s",iteend,chout);
  system(commande);
  sprintf(commande,"mv res.%s REJECTED/%s",chout,chout);
  system(commande);         
  sprintf(commande,"mv DEPHA* REJECTED/%s",chout);
  system(commande);         
  sprintf(commande,"mv energiemode* REJECTED/%s",chout);
  system(commande); 
  sprintf(commande,"mv error? ??.????.??.lhz REJECTED/%s",chout);  
  system(commande);  
  sprintf(commande,"mv AMPL_RATIO REJECTED/%s",chout);  
  system(commande);  
}  
sprintf(commande,"gzip %s %s %s ",yiray,souray,dpray); 
system(commande);   
sprintf(commande,"gzip %s %s %s ",yilov,soulov,dplov); 
system(commande);   
/*system ("rm error? *.it* des.so.? des.so.??");*/
system ("rm *.it* des.so.? des.so.??");
system("rm *.init*");
system("rm *.postscript");
/* system("rm initrun newcorp autonewrun bid ECH"); */
system("rm initrun newcorp autonewrun bid" ); 
system("rm DEPH* energiemode* pick* res*");
system("rm *.lhz");
}
/*###################################################################################*/
void initnewrun(istart)
int istart;
{
FILE * out;

system ("rm bid");
system ("rm newcorp");

if(!(out=fopen("bid","w")))
    {printf ("\n pb opening file bid");
     exit(1);}
fprintf(out,"@ itstart  =  %d",istart);
fprintf(out,"\nset FicTV  =  tv.it%d",istart);
fprintf(out,"\nset DISPray =     disp.ray.it%d",istart);
fprintf(out,"\nset DISPlov =     disp.lov.it%d",istart);
fprintf(out,"\nset ATTray  =     att.ray.it%d",istart);
fprintf(out,"\nset ATTlov  =     att.lov.it%d",istart);
fprintf(out,"\nset GKray1  =     gk.ray.1.it%d",istart);
fprintf(out,"\nset GKlov1  =     gk.lov.1.it%d",istart);
fprintf(out,"\nset TFray1  =     tf.ray.1.it%d",istart);
fprintf(out,"\nset TFlov1  =     tf.lov.1.it%d\n",istart);
fclose(out);
system ("cat bid corpnewrun >newcorp");
system ("rm bid");
}
/*######################################################################################*/
void initperiode(dista,periode,idata,ilov)
float dista;
int ilov;
int *periode;       
int idata;           
{
int compt,compt1,nper,npertot,i,ii;
int per[4],pertot[4];
char filenm[20];
int f40,nmode;
FILE *tper;

/* cette routine choisit les periodes de filtrage au
   vu des info lues dans les fichiers DATATOTRAY et DATATOTLOV
   les ecrit dans un fichier pick.ray.1 ou pick.lov.1 qui seront lu
   dans cross12.lr.auto. Le choix des periodes est base sur le fait que
   dans preselec.t.c on ne selectionne un sismogramme que si il presente
   un bon rapport signal/bruit au moins a 30 et 60s ou 40 et 80s.
   Ici la selection des periodes est faite selon l'ordre de preference suivant :
   40-80-160  1er choix
   30-60-120  2eme choix
   40-80      3eme choix 
   30-60      4eme choix.
   Les possibilites 80-160, 60-120, 30, 40, 60, 80, 120, 160 doivent theoriquement
   n'etres jamais rencontre. Si c'est le cas on arrete la routine et passe au
   seisme suivant
   --------------------------------------------------------------------
   Modif 24 Octobre 1997 suite une suggestion de Brian
   30s c'est trop courte periode pour du fondamental en continental.
   Au lieu de piquer les enveloppes a 30s on les piques a 40
   les choix deviennent :
   40-80-160  1er choix
   40-60-120  2eme choix
   40-80      3eme choix
   40-60      4eme choix
   --------------------------------------------------------------------
   Modif 8 Decembre 1997 : On se donne la possibilite de ne pas inverser
   les periodes >100s (donnees SKIPPY ou deuxieme run des LP).
   --------------------------------------------------------------------
   Modif 25 Mars 2000 : On se donne la possibilite de travailler a plus 
   longue periode 
   50-90-160  1er choix
   50-70-120  2eme choix
   50-90      3eme choix
   50-70      4eme choix
    
 */

nmode=5;
tmin=0;
tmax=0;
compt=0;
/* Modif 10 Juillet 1998 : on prevoit le cas ou les criteres de  
  rapport signal/bruit pour retenir un sismo dans preselec.t sont
  effectues a 40-60-120 au lieu de 30-60-120 : c'est ce qui a ete fait
  lors de la deuxieme selection des donnees SKIPPY, etant donne qu'on avait
  fait le choix de ne plus inverser le 30s.  */
/*if (periode[0]==1)compt=1;*/
if ( (periode[0]==1) || (periode[1]==1) )compt=1;
if (periode[2]==1)compt=compt+3;
if (periode[4]==1)compt=compt+5;

compt1=0;
if (periode[1]==1)compt1=1;
if (periode[3]==1)compt1=compt1+3;
if (periode[5]==1)compt1=compt1+5;

if  (((compt==9) || (compt==4)) &&
     ((compt1==9) || (compt1==4)))
    {
     if(compt1>=compt)
         { compt=compt1;
           f40=1;}
     else if(compt1<compt)
         { f40=0;}
    } 
else if (((compt==9) || (compt==4)) &&
         ((compt1!=9) && (compt1!=4)))
    {
     f40=0;
    }
else if (((compt!=9) && (compt!=4)) &&    
         ((compt1==9) || (compt1==4)))
    { 
     compt=compt1;
     f40=1;
    }

if (((compt==1) || (compt==3) || (compt==5) || (compt==6) || (compt==8)) &&
    ((compt1==1) || (compt1==3) || (compt1==5) || (compt1==6) || (compt1==8)))
     {printf ("\n Ce pgm suppose que le sismogramme est modelise seult si:");
      printf ("\n -f40 et f80 sont ok ou f30 et f60 sont ok (pgm preselec.t.c)");
      printf ("\n ou (modif 10 Juillet 98 si f40 et f60 ok )");
      printf ("\n Ce n'est visiblement pas le cas ici=>beug quelque part.");
      printf("\n On a compt = %d and compt1 = %d",compt,compt1);
      exit(1);}

pertot[0]=0; 
pertot[1]=0; 
pertot[2]=0; 
pertot[3]=0; 

/* Modif 8 Decembre 1997 : on inverse les periodes >100s au premier run. 
   On fera tourner l'algorithm une deuxieme
   fois sur les donnees rejetees en inversant uniquement donnees a T<100s. */ 

/* if (compt==9) compt=4; */
 if (compt==4) 
    {    
     npertot=2;   
     if (f40==0) 
          {pertot[0]=50;  /* pertot[0]=30; */
           pertot[1]=70;} 
     if (f40==1) 
          {pertot[0]=50;  
           pertot[1]=90;} 
     if ((f40==0) && (idata==-2))
          {per[0]=70; 
           nper=1;   }
     else if ((f40==0) && (idata>=-1)) 
          {per[0]=50;  /*per[0]=30; */
           per[1]=70; 
           nper=2;}
     else if ((f40==1) && (idata==-2)) 
          {per[0]=90; 
           nper=1;}
     else if ((f40==1) && (idata>=-1)) 
          {per[0]=50; 
           per[1]=90;
           nper=2;}
    }

if (compt==9)
    {
     npertot=3;   
     if (f40==0) 
          {pertot[0]=50;  /*pertot[0]=30; */
           pertot[1]=70;  
           pertot[2]=120;} 
     if (f40==1) 
          {pertot[0]=50;  
           pertot[1]=90; 
           pertot[2]=160;} 
     if ((f40==0)  && (idata==-2))
            {per[0]=70; 
             per[1]=120; 
             nper=2;}
     else if ((f40==0)  && (idata>=-1))
            {per[0]=50;  
             per[1]=70; 
             per[2]=120;
             nper=3;}
     else if ((f40==1) && (idata==-2))
            {per[0]=90; 
             per[1]=160; 
             nper=2;}
     else if ((f40==1) && (idata>=-1))
            {per[0]=50; 
             per[1]=90; 
             per[2]=160;
             nper=3;}
    } 

if(ilov==0) strcpy(filenm,"pick.ray.1");
else if (ilov==1) strcpy(filenm,"pick.lov.1");

/* 
on considere que l'onde de surface est formee a une 
periode qd dist >=10*periode                    
*/

printf ("\n NPER=%d",nper);
for (ii=0;ii<nper;ii++)          
    if(dista<=(float)(10*per[ii])) nper-- ;

if (nper == 0)
  {
   printf ("\n Pas de periodes utilisables pour ce trajet????");
   printf ("\n D'apres les criteres de seections, ie :");
   printf ("\n distance >1000km et 2 T, 30-60 ou 40 80 disponibles");
   printf ("\n ce cas ne devrais jamais etre rencontre\n");
   printf ("\n DIST=%f",dista);
   exit(1);
  }

if (f40==0) tmin=50; /*tmin=30;*/
else if (f40==1) tmin=50;
tmax=per[nper-1];

if(!(tper=fopen(filenm,"w")))
  {printf ("\n pb opening file : %s",filenm);
   exit(1);}
fprintf(tper," %d %d %d %d %d",nmode,npertot,pertot[0],pertot[1],pertot[2]);
fprintf(tper,"\no");
  for (ii=0;ii<nper;ii++)
   {
    fprintf(tper,"\no");
    fprintf(tper, "\n %d",per[ii]);
    printf( "\n PER %d %d",ii,per[ii]);
    fprintf(tper,"\n 1");
   
/*modif le 05/06/2002 : le dernier chiffre ecrit dans pick.* doit etre un 
  interger puisque cross12 essaie ensuite de lire un interger. Pas de pb
  avec unix mais bug avec linux */

    if(((ii+1)==nper) && (idata==-2)) 
       {fprintf(tper,"\nn");
        fprintf(tper,"\n 0");}
    else if(((ii+1)==nper) && (idata==-1))
       {fprintf(tper,"\nn");
        fprintf(tper,"\n 0");}
    else if(((ii+1)==nper) && (idata==0))
       {fprintf(tper,"\no");
        fprintf(tper, "\n %d",per[1]);
        fprintf(tper,"\n 1");
        fprintf(tper,"\nn");
        fprintf(tper,"\n 1");}
    else if(((ii+1)==nper) && (idata==1))
       {fprintf(tper,"\no");
        fprintf(tper, "\n %d",per[1]);
        fprintf(tper,"\n 1");
        fprintf(tper,"\no");
        fprintf(tper, "\n %d",per[0]);
        fprintf(tper,"\n 1");
        fprintf(tper,"\nn");
        fprintf(tper,"\n 2");}
    }
    fprintf(tper,"\n");
fclose(tper);
}

/*#################################################################################*/
void initialisation(path_data,path_cmt,path_saito,path_resp,nomdatray,nomdatlov,adrdatray,adrdatlov,ilr,yiray,yilov,souray,soulov,dpray,dplov)
char adrdatray[80],adrdatlov[80];
char path_data[150],path_cmt[150],path_saito[150],path_resp[150];
char nomdatray[45],nomdatlov[45];
char * yiray,* yilov,* souray,* soulov,* dpray,* dplov;
int ilr;
/* 13 octobre 98 : legere modif pour inverser Love/Rayleigh simultanement.
   J'ai cree un adrsaitoray et un adrsaitolov qui permettent de stocker les
   derivees partielles et deplacemnet tension Love dans des directory differents de
   ceux utilises pour Rayleigh
   9 Aout 99 : j'ai rajoute un adrsaitosouray et un adrsaitosoulov qui permet de
   stocker les fichiers sources dans des directorys differents. Pour travailler avec
   prem comme modele initial, l'extension "prem" est ajoutee a l'aide d'un strcpy*/
{
struct cmt *cmt;
FILE *fichdat;
FILE *sortie;
float toto;
int i,nbsei,nlr,seiray,npts;
int yearo,montho,dayo,heuro,mino,bidi;
int tzmaxray, tzminray;
int tzmaxlov, tzminlov;
float seco,bidf,distrad;
float elat,elon,evdp,slat,slon,salt,pas;
char ligne[100],instray[20],instlov[20],adrsaitoray[150],adrsaitolov[150],adrsaitosouray[150],adrsaitosoulov[150],bids[20];
char commande[300], foo[10],resp[20];

if ((cmt=(struct cmt*)malloc (sizeof(struct cmt)))==NULL)
     {printf ("\n Error in allocation struct cmt");
      exit(1);}
/* ouverture du fichier dat rayleigh -> lat,lon,dist,azim */
   sprintf(adrdatray,"%s%s",path_data,nomdatray);      
/*   printf("nomdatray, adrdatray is %s, %s\n", nomdatray, adrdatray); */
   
if(!(fichdat=fopen(adrdatray,"r")))
       {printf ("\n cannot open file bar : %s\n",adrdatray);
        exit(1);} 
fgets(ligne,100,fichdat);
fgets(ligne,100,fichdat);

fgets(ligne,100,fichdat);
sscanf(ligne,"%s %s %f %s %f %s %f %s",bids,bids,&elat,bids,&elon,bids,&toto,bids);
cmt->evdp=toto;

fgets(ligne,100,fichdat);
sscanf(ligne,"%s %s %f %s %f %s %f %s",bids,bids,&slat,bids,&slon,bids,&salt,bids);

fgets(ligne,100,fichdat);
sscanf(ligne,"%d %d %f",&heuro,&mino,&seco);

fgets(ligne,100,fichdat);
sscanf(ligne,"%d %d %f %s %d %d %d %s",&bidi,&bidi,&bidf,bids,&dayo,&montho,&yearo,bids);

fgets(ligne,100,fichdat);
sscanf(ligne,"%d %f %f %s",&npts,&pas,&toto,bids);
cmt->dist=toto;

/* printf("\n%f %f %f %f %f %f\n",elat,elon,slat,slon,cmt->evdp,cmt->dist); */
fclose (fichdat);
elat=atan(0.993277*tan(elat*pi/180.));
slat=atan(0.993277*tan(slat*pi/180.));
elon=elon*pi/180.;
slon=slon*pi/180.;
distrad=(cmt->dist/111.195)*pi/180.;
cmt->azim=azimuth(elat,elon,slat,slon,distrad);
cmt->azim=cmt->azim*180./pi;
/* printf("\n CMT-AZIM %f ",cmt->azim); */

/* recuperation dans le fichier CMT des parametres a la source*/
readcmt(path_cmt,nomdatray,cmt,yearo,montho,dayo,heuro,mino,seco);
   printf("hello4?");

   strcpy (adrsaitosoulov,path_saito);
   strcpy (adrsaitosouray,path_saito);
   strcpy (adrsaitolov,path_saito);
   strcpy (adrsaitoray,path_saito);

   printf("hello?");
	   
/* strcpy (adrsaito,"./"); */

 /* If string IU.PABB....LHZ is contained in the data file name, then set response to IU.PABB....lhz. This way you mush the data file to the right response */  
/* MUST PUT HERE MY STATIONS  */
/* WHY not cutting the filename and create the IU....lhz name???? -- MUCH MORE GENERAL*/


/* Sylvana's change:  Get response file name from *.dat file name instead of specifying it explicitly in the code */

 strcpy(resp,"");

  for (i=23; i<=33; i=i+1) {
    sprintf(foo, "%c",nomdatray[i]);
    strcat(resp,foo);
    printf("RESP IS %s\n",resp);
  }

	  
  /* add lhz at end of name */
  /* rayleigh */
  strcpy (instray,resp);
  strcat(instray,"lhz");
     /* love */
  strcpy (instlov,resp);
  strcat(instlov,"lhn");
	  
  printf("response names are %s %s\n", instray,instlov);
  
/*if(!(strstr(nomdatray,"IU.PABB....LHZ")==NULL))
    {strcpy(instray,"IU.PABB....lhz");
     strcpy(instlov,"IU.PABB....lhn");}
else
   {printf ("\n Quel instrument pour %s ",nomdatray);
    exit(1);} */

/* ------------------------------------------------- */

strcpy(dpray,adrsaitoray);
strcat(dpray,"dp.mod.");
strcat(dpray,nomdatray);
strcat(dpray,".SPR");
printf(" FILE DPRAY IS %s\n", dpray);

strcpy(yiray,adrsaitoray);
strcat(yiray,"yi.mod.");
strcat(yiray,nomdatray);
strcat(yiray,".SPR");
printf(" FILE YIRAY IS %s\n", yiray);

strcpy(souray,adrsaitosouray);
strcat(souray,"yi.mod.sour.");
strcat(souray,nomdatray);
printf(" FILE SOURAY IS %s\n", souray);

sprintf(commande,"gunzip %s %s %s",yiray,souray,dpray);
system(commande);

printf("2\n");  
printf(" FILE DPRAY IS %s\n", dpray);
printf(" FILE YIRAY IS %s\n", yiray);
printf(" FILE SOURAY IS %s\n", souray);

  if(tpminray>=5)tzminray=tpminray-5;  /* tpmin est normalement 30 ou 40: on regarde le sismo aux */
  else tzminray=0;
  if(tpmaxray>=5)tzmaxray=tpmaxray*2;  /* periodes ou on l'inverse */
  else tzmaxray=0;
  if(tpminlov>=5)tzminlov=tpminlov-5;
  else tzminlov=0;
  if(tpmaxlov>=5)tzmaxlov=tpmaxlov*2;
  else tzmaxlov=0;
/* affectation des autres parametres*/
printf("3\n");  
printf(" FILE DPRAY IS %s\n", dpray);
printf(" FILE YIRAY IS %s\n", yiray);
printf(" FILE SOURAY IS %s\n", souray);


if(ilr==0)
{
nbsei=1;
nlr=0;
seiray=1;
strcpy(dplov,"bid");
strcpy(yilov,"bid");
strcpy(soulov,"bid");
}
else if(ilr==1)
{
nbsei=1;
nlr=1;
seiray=1;
strcpy(dplov,adrsaitolov);
strcat(dplov,"dp.mod.");
strcat(dplov,nomdatlov);
strcat(dplov,".SPR");
strcpy(yilov,adrsaitolov);
strcat(yilov,"yi.mod.");
strcat(yilov,nomdatlov);
strcat(yilov,".SPR");
strcpy(soulov,adrsaitosoulov);
strcat(soulov,"yi.mod.sour.");
strcat(soulov,nomdatlov);
sprintf(commande,"gunzip %s %s %s",yilov,soulov,dplov);
system(commande);
}


printf("lov\n");  
printf(" FILE DP IS %s\n", dplov);
printf(" FILE YI IS %s\n", yilov);
printf(" FILE SOUR IS %s\n", soulov);

printf("4\n");  
printf(" FILE DPRAY IS %s\n", dpray);
printf(" FILE YIRAY IS %s\n", yiray);
printf(" FILE SOURAY IS %s\n", souray);

/*ecriture des affectations dans le fichier de sortie*/
if(!(sortie=fopen("initrun","w")))
       {printf ("\n pb opening file %s ",nomdatray);
        exit(1);}
fprintf(sortie,"#!/bin/csh");
fprintf(sortie,"\n@ NBSEI = %d", nbsei);
fprintf(sortie,"\n@ NLR = %d",nlr);
fprintf(sortie,"\n@ SEIray =  %d",seiray);
fprintf(sortie,"\nset DONray1   =     %s",adrdatray);
fprintf(sortie,"\nset DONlov1   =     %s",adrdatlov);
fprintf(sortie,"\nset DERPARTray = %s",dpray);
fprintf(sortie,"\nset DERPARTlov = %s",dplov);
fprintf(sortie,"\nset DEPTENSray =  %s",yiray);
fprintf(sortie,"\nset DEPTENSlov =  %s",yilov);
fprintf(sortie,"\nset DEPTSOUray = %s",souray);
fprintf(sortie,"\nset DEPTSOUlov = %s",soulov);
fprintf(sortie,"\nset INSTray = %s%s",path_resp,instray);
fprintf(sortie,"\nset INSTray2 = %s",instray);
fprintf(sortie,"\nset INSTlov = %s%s",path_resp,instlov);
fprintf(sortie,"\nset INSTlov2 = %s",instlov);
fprintf(sortie,"\nset AXEP1     =  '%d %d'",cmt->axepp,cmt->axepa);
fprintf(sortie,"\nset AXET1     =  '%d %d'",cmt->axetp,cmt->axeta);
fprintf(sortie,"\nset PROF1     =   %f",cmt->evdp);
fprintf(sortie,"\nset MOMENT1   =   %E",cmt->moment);
fprintf(sortie,"\nset DIST1     =   %f",cmt->dist);
fprintf(sortie,"\nset CENTROID1 =  oui");
fprintf(sortie,"\nset AZIM1     =   %f",cmt->azim);
fprintf(sortie,"\nset DUREE1    =   %f",cmt->dur);
fprintf(sortie,"\nset RISETIME1 =  0.");
fprintf(sortie,"\nset TPLATEAUray1 = '%d %d'",tpminray,tpmaxray);
fprintf(sortie,"\nset TRETOUR0ray1 = '%d %d'",tzminray,tzmaxray);
fprintf(sortie,"\nset TPLATEAUlov1 = '%d %d'",tpminlov,tpmaxlov);
fprintf(sortie,"\nset TRETOUR0lov1 = '%d %d'",tzminlov,tzmaxlov);
fprintf(sortie,"\nset two = '2'");
fprintf(sortie,"\n");
fclose(sortie);
}
/*#################################################################################*/
void readcmt(path_cmt,nomdatray,cmt,yearo,montho,dayo,heuro,mino,seco)
char nomdatray[45];
char path_cmt[150];
int yearo,montho,dayo,heuro,mino;
float seco;
struct cmt *cmt;
{
FILE *entree;
char ligne1[200],ligne2[200],ligne3[200],ligne4[200];
char bidc;
int ico;
/* cette routine est censee renvoyer a l'initialisation les parametres*/
/* a la source lus dans le fichier CMT-SELECTED qui sera copie*/
/* dans le directory d'utilisation pour les tests histoires d'eviter les*/
/* grosses conneries.*/
/* Declaration de variables dans le fichier CMT-SELECTED  */
/*-----ligne 1--------------------------------------------*/
       char code[9],region[25],sousligne[50];
       int month,day,year,hour,min;
       float sec;
       float lat,lon,depth,mb,ms;
       int val1,pd1,val2,pd2,val3,pd3;
/*-----ligne 2--------------------------------------------*/
       char epsrc[8];
       float edt,eclat,eclon,ecdepth;
       float dt,clat,clon,cdepth;
       int bwst,bwrec,bwcoff;
       int mwst,mwrec,mwcoff;
/*-----ligne 3--------------------------------------------*/
       float ahdur,mrr,emrr,mss,emss,mee,emee;
       float mrs,emrs,mre,emre,mse,emse;
       int expo;
       char bids[8];
/*-----ligne 4--------------------------------------------*/
       float smo, ev[4];
       int evpl[4],evaz[4],strike[3],dip[3],rake[3];


if(!(entree=fopen(path_cmt,"r")))
/* if(!(entree=fopen("/home/pilidou/programs/codes_CMT_sorter_Ms/test","r"))) */
       {printf ("\n cannot open file %s", path_cmt);
        exit(1);}

  while(!feof(entree))
  {
    fgets(ligne1,200,entree);
    if(feof(entree)) goto sortieboucle;
    fgets(ligne2,200,entree);
    fgets(ligne3,200,entree);
    fgets(ligne4,200,entree);

/*  sscanf(ligne1,"%9s %2d %1s %2d %1s %2d %3d %1s %2d %c %f %f %f %3d %c %1d %1d %c %1d %1d %c %1d %s"*/
    sscanf(ligne1,"%9s %2d%c%2d%c%2d %3d%c%2d%c%f %f %f %3d%c%1d%1d%c%1d%1d%c%1d%s"
   ,code,&month,bids,&day,bids,&year,&hour,bids,&min,&bidc,&sec,&lat,&lon,&val1
   ,&bidc,&pd1,&val2,&bidc,&pd2,&val3,&bidc,&pd3,region);
    depth=(float)val1+(float)pd1/10.;
    mb=(float)val2+(float)pd2/10.;
    ms=(float)val3+(float)pd3/10.;

/*  sscanf(ligne2,"%s%3s%2d%3d%4d%3s%2d%3d%4d %s %f %f %f %f %f %f %f %f" */
    sscanf(ligne2,"%s%c%c%c%c%2d%3d%4d%c%c%c%c%2d%3d%4d%c%c%c%c %f %f %f %f %f %f %f %f"
   ,epsrc,&bidc,&bidc,&bidc,&bidc,&bwst,&bwrec,&bwcoff,&bidc,&bidc,&bidc,&bidc
   ,&mwst,&mwrec,&mwcoff,&bidc,&bidc,&bidc,&bidc,&dt,&edt,&clat,&eclat
   ,&clon,&eclon,&cdepth,&ecdepth);

    sscanf(ligne3,"%c%c%c%c %f%c%c%c%c%2d %f %f %f %f %f %f %f %f %f %f %f %f"
    ,&bidc,&bidc,&bidc,&bidc,&ahdur,&bidc,&bidc,&bidc,&bidc,&expo
    ,&mrr,&emrr,&mss,&emss,&mee,&emee,&mrs,&emrs
    ,&mre,&emre,&mse,&emse);
    
    sscanf(ligne4," %f %d %d %f %d %d %f %d %d %f %d %d %d %d %d %d"
     ,&ev[1],&evpl[1],&evaz[1]
     ,&ev[2],&evpl[2],&evaz[2]
     ,&ev[3],&evpl[3],&evaz[3],&smo
     ,&strike[1],&dip[1],&rake[1]
     ,&strike[2],&dip[2],&rake[2]);
/*  
    printf("\n\n %s %d %d %d %d %d %f %f %f %f %f %f %s"
   ,code,month,day,year,hour,min,sec,lat,lon,depth,mb,ms,region);
    printf("\n %s %d %d %d %d %d %d %f %f %f %f %f %f %f %f" 
   ,epsrc,bwst,bwrec,bwcoff,mwst,mwrec,mwcoff,dt,edt,clat,eclat, 
    clon,eclon,cdepth,ecdepth); 
    printf("\n %f %d %f %f %f %f %f %f %f %f %f %f %f %f"
    ,ahdur,expo,mrr,emrr,mss,emss,mee,emee,mrs,emrs,mre,emre,mse,emse);
    printf("\n %f %d %d \n%f %d %d \n%f %d %d %f \n%d %d %d \n%d %d %d"
     ,ev[1],evpl[1],evaz[1]
     ,ev[2],evpl[2],evaz[2]
     ,ev[3],evpl[3],evaz[3],smo
     ,strike[1],dip[1],rake[1]
     ,strike[2],dip[2],rake[2]);
*/
    
/*
1001 format(a8,5(1x,i2),1x,f4.1,f7.2,f8.2,f6.1,2f3.1,a24)
1002 format(a3,2(4x,i2,i3,i4),4x,f6.1,f4.1,
     1       f7.2,f5.2,f8.2,f5.2,f6.1,f5.1)
1003 format(4x,f4.1,4x,i2,6(f6.2,f5.2))
1004 format(3(f7.2,i3,i4),f7.2,2(i4,i3,i5))
*/

  /* Remise a jour du temps origine lu dans CMT-SELECTED avec
     celui present dans le fichier .dat qui est le centroid (tempsCMT* +DT)*/

sec=sec+dt;
if(sec>=60.)
{min++;
 sec=sec-60.;}
else if (sec<0.)
{min--;
 sec=sec+60.;}

if(min>=60)
{hour++;
 min=min-60.;}
else if (min<0.)
{hour--;
 min=min+60.;}

if(hour>=24)
{day++;
 hour=hour-24;}
else if (hour<0.)
{day--;
 hour=hour+24;}

if((day>365) || (day<0)) 
 {printf ("Peculiar Earthquake not compatible with sortdata1bis\n");
  exit(1);}
  /*
    printf("\n ");
    printf("\n %d %d %d %d %d %f",year,month,day,hour,min,sec);
    printf("\n ");
    */

/* --------------------------------------------------------------- */
  /* Sylvana change: added years upto 2020 !!!!, only accounted for 2000 and code was making year 2001 1901 etc ==> silly hige numbers for any cmt allorder values. */ 
if (year < 20) 	year=year+2000;
else year=year+1900;
/* --------------------------------------------------------------- */
/*
    printf("\n %d %d %d %d %d %f",yearo,montho,dayo,heuro,mino,seco);
    printf("\n %d %d %d %d %d %f",year,month,day,hour,min,sec);
*/    
if((year==yearo) && (month==montho) && (day==dayo) && (hour==heuro) && (min==mino) && (fabs(sec-seco)<1.0))
      {
	printf("\n %d %d %d %d %d %f",yearo,montho,dayo,heuro,mino,seco);
	printf("\n %d %d %d %d %d %f",year,month,day,hour,min,sec);
        cmt->axetp=evpl[1];
        cmt->axeta=evaz[1];
        cmt->axepp=evpl[3];
        cmt->axepa=evaz[3];
        cmt->moment=smo*pow(10.,(float)expo); /* verifie : cf fich CMT.FORMATEXPLANATION */
        cmt->dur=ahdur;
      }
 }
sortieboucle:;
printf("\n plongement azimuth axet %d  %d",cmt->axetp,cmt->axeta);
printf("\n plongement azimuth axep %d  %d",cmt->axepp,cmt->axepa);
printf("\n Moment sismique %E",cmt->moment);
printf("\n 1/2 duration %f",cmt->dur);
printf("\n FIN LECTURE ");
printf("\n hello0?");
/*exit(1); */
fclose(entree);
printf("\n hello1?");
}
/* #########################great cicle azimuth####################################### */
/*The argument late,lone,lats,lons,delt are supposed to be given
  in radians. late and lats are latitude (not colatitude)*/
float azimuth(late,lone,lats,lons,delt)
float late,lone,lats,lons,delt;
{
#define Epsilon 0.00000001
#define pi 3.14159265358979323846264338
double Cae,Sae,Cas,Sas,Clt,Slt,Sll;
double Caz,Saz, azim;    

Cae=cos(late);   
Sae=sin(late);
Cas=cos(lats);
Sas=sin(lats);
Sll=sin(lons-lone);
Clt=cos(delt);    
Slt=sin(delt);

Caz=(Sas-Sae*Clt)/(Cae*Slt);
Saz=(Cas*Sll)/Slt;

if((Saz/Caz)>=0.) azim=atan(Saz/Caz);
else azim=atan(-(Saz/Caz));

/*printf ("\n azim : %f",azim);*/

        if (Saz>=0.)
         {
          if (Caz>0)
            azim=azim;
          else
            azim=(pi-azim);
         }   
        else
         { 
          if (Caz>0.)
            azim=(2*pi-azim);
          else
            azim=(pi+azim);
         }   
azim=(float)azim;
return (azim);
}
/*#############################################################################*/
