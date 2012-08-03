#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
/*######Cette version ecrit la liste des fichiers traites dans lst.mod####*/
/* Code by Eric Debayle*/
/* 23 Mai 2002 : Le code a ete modifie pour :
     - prendre en compte l'epaisseur de glace
     - ecrire les fichiers mod dans un format compatible avec saito02.f
       (2 couches d'epaisseur non nulles entre chaque discontinuite dans la croute)
   la structure crustale comprend maintenant une couche d'eau, une de glace, une de sediment,
   une de croute sup et une de croute inf. Si les couches d'eau, de glace ou de sediment sont nulles
   elles sont remplacees par une couche de croute sup pour garder un nombre de couche constant.
*/
/* 15 juin 2000 : tout les problemes d'arrondis lors de l'ecriture des epaisseurs
   (au moins tout ceux detectes ont ete evites)
   21 Mars 2000 : cette version cree des fichiers mod avec la croute 3SMAC
   et le modeles prem lisse Qcst pour le manteau (dans /home/eric/3-SMAC.chantier/NEWLR) */
/* modifie le 23 juillet 97=> longitude ramenee entre 0 et 360.
   on prend un Q de 600 pour la croute qqsoit la valeur lue dans 3SMAC
   les valeurs de 3SMAC apparaissant particulierement faible a la base
   de la croute sous l'Australie */
/* ##### Quelques rappels 3SMAC ########### */
/*   3SMAC : une donnee au centre de carres de 2*2degres
     180 points en longitudes en partant de 0degres
     90 points en lat en partant du pole Nord.
     Les indices i,j ont ete verifies , sont ok.    */
/* #####  Facteur de qualite dans 3SMAC ########### */
/*   Qmu: c'est la meme quantite physique  que Qbeta
     on a:   Qmu = Qbeta */
/* ##### VARIABLES GLOBALES ########### */
float lati,loni;
#define pi 3.14159265358979323846264338

/* ##### declaration de structure ########*/
struct layer
{
/*NB : 400 est choisi pour maximiser le nombre de cellules
       de 2*2 degres echantillonables le long d'un trajet */
   
float epe[400];
float epg[400];
float eps[400];
float epcs[400];
float epci[400];
float vp[10][400];
float vs[10][400];
float rho[10][400];
float Q[10][400];
};
struct resu 
{
float ep[10];
float vp[10];
float vs[10];
float rho[10];
float Q[10];
};
struct source
{
float ep[10];
float vp[25];
float vs[25];
float rho[25];
float Q[25];
};
struct smac
{
float ep[5][180][90];
float vp[10][180][90];
float vs[10][180][90];
float rho[10][180][90];
float Qmu[10][180][90];
float ice[180][90];
};
struct mant
{
float vp[15][180][90];
float vs[15][180][90];
float rho[15][180][90];
float Qmu[15][180][90];
};
/*#####routine declaration#######*/
void moyo();
void affstructsmac();
float ** lecsmac();
float ** allotab();
float del();
float azimuth();
void coord();
void woutave();
void woutsou();
/* ####  MAIN PGM   ################### */
main()
{
float elat,elon,slat,slon,azes;
float dpas,delta,dstcle;
float coelat,coslat;
int npts[400];
int nt,il,i,j,ii,jj,ic,imoins,jmoins,indbloc,ipath,lr,npaths;
struct smac *smac;
struct mant *mant;
struct layer *tra;
struct resu *resu;
struct source *sour;
char nomfich[100];
char nomsort[100];
char nomso[100];
char nombid[100];
char nomin[100];
char nomin2[100];
char ligne[100];
FILE *lssortie;

  /*memory allocation for structure smac and mant*/
  if ((smac=(struct smac*)malloc (sizeof(struct smac)))==NULL)
     {printf ("\n Error in allocation struct smac");
      exit(1);}
  if ((mant=(struct mant*)malloc (sizeof(struct mant)))==NULL)
     {printf ("\n Error in allocation struct mant");
      exit(1);}
  affstructsmac(smac, mant);

  printf ("\n number of paths for which you want a crustal model?\t");
  gets(ligne);
  sscanf(ligne,"%d ",&npaths);

  strcpy(nomsort,"lst.mod");
  if(!(lssortie=fopen(nomsort,"w")))
       {printf ("%s\n : pb opening file",nomsort);
        exit(1);}
  fprintf (lssortie, "#");
  fprintf (lssortie, "\n set NBSEI = %d", npaths);

  for (ipath=1;ipath<=npaths;ipath++)
  {
      printf ("\n name output file, epicenter and station coord(lat,lon) ?\t");
      gets(ligne);
      sscanf(ligne,"%s %f %f %f %f ",nomsort,&elat,&elon,&slat,&slon);
      if(elon<0) elon=elon+360.;
      if((elon<0.) || (elon>360.)) 
        {printf ("\n Pb with lon path %s",nomsort);
         exit(1);}
    
  /*memory allocation for structure layer and sour*/
  if ((tra=(struct layer*)malloc (sizeof(struct layer)))==NULL)
     {printf ("\n Error in allocation struct layer");
      exit(1);}
     
  if ((sour=(struct source*)malloc (sizeof(struct source)))==NULL)
     {printf ("\n Error in allocation struct layer");
      exit(1);}
     
      elat=atan(0.993277*tan(elat*pi/180.));
      slat=atan(0.993277*tan(slat*pi/180.));
      coelat=(pi/2)-elat;
      coslat=(pi/2)-slat;
      elon=elon*pi/180;
      slon=slon*pi/180;
      delta=del(coelat,elon,coslat,slon);
      azes=azimuth(elat,elon,slat,slon,delta);
      azes=azes*180/pi;
      delta=delta*180/pi;
    
/*    printf ("\n azimuth et distance en degres");
      printf ("\n %f %f",azes,delta); */
    
      if (azes<0  || azes>360)
         {printf ("\n azimuth<0 ou>360");
          exit(1);}
      nt=(int)(5.*delta+1.);
      dpas=delta/(float)nt;
      nt++;

    /* initialisation   */                     
    
      for(jj=0;jj<400;jj++) npts[jj]=0.;

      i=9999999;
      j=9999999;
      indbloc=0;
      
    /* attention les indices i et j vont de 0 a 179
       et de 0 a 89 */
    
     printf("\n NT %d",nt);  
      for(jj=1;jj<=nt;jj++)
        {
         dstcle=dpas*(float)(jj-1);
         coord(coelat,elon,dstcle*pi/180,azes*pi/180);
     /*  printf ("\n LATI %f LONI %f ",lati,loni);  */
         imoins=i;
         jmoins=j;
         if(loni<0.) loni=loni+360.;
         if(loni>360.) loni=loni-360.;
         if((loni<0.) || (loni>360.)) 
           {printf ("\n Pb with computed long along path %s",nomsort);
            exit(1);}
         i=(int)(loni/2.); 
         j=(int)((90.-lati)/2.);

         /*on recupere le modele crustal et mantellique 
           sous l'epicentre */
         if (jj==1)
            {
            for(ii=0;ii<5;ii++) 
             { sour->ep[2*ii+1]= smac->ep[ii][i][j];
               sour->ep[2*ii+1]= (int)(100*sour->ep[2*ii+1])/100.; } /* <- avoid pbs when printing ep*/
            for(ii=0;ii<10;ii++) 
               {
               sour->vp[ii]= smac->vp[ii][i][j];
               sour->vs[ii]= smac->vs[ii][i][j];
               sour->rho[ii]= smac->rho[ii][i][j];
               sour->Q[ii]= smac->Qmu[ii][i][j];
               }
            for(ii=0;ii<15;ii++) 
               { 
                ic=ii+10; 
                sour->vp[ic]= mant->vp[ii][i][j];  
                sour->vs[ic]= mant->vs[ii][i][j];
                sour->rho[ic]= mant->rho[ii][i][j];
                sour->Q[ic]= mant->Qmu[ii][i][j];
               }   
            }  

         if ((i==imoins)&&(j==jmoins))
            {
            npts[indbloc]++;
            }
         else
            {
               indbloc++;
               npts[indbloc]=1;
               tra->epe[indbloc]=smac->ep[0][i][j];
               tra->epg[indbloc]=smac->ep[1][i][j];
               tra->eps[indbloc]=smac->ep[2][i][j];
               tra->epcs[indbloc]=smac->ep[3][i][j];
               tra->epci[indbloc]=smac->ep[4][i][j];
       /*    printf ("\n bloc %d  %f %f %f %f %f %f %f",indbloc,lati,loni,tra->epe[indbloc],tra->epg[indbloc],tra->eps[indbloc],tra->epcs[indbloc],tra->epci[indbloc]);  */
            for(ii=0;ii<10;ii++)
               {
               tra->vp[ii][indbloc]=smac->vp[ii][i][j];
               tra->vs[ii][indbloc]=smac->vs[ii][i][j];
               tra->rho[ii][indbloc]=smac->rho[ii][i][j];
               tra->Q[ii][indbloc]=smac->Qmu[ii][i][j];
               }
             }
    /*     printf ("\n I %d I %d ",i,j);  */
          }
 /*  for(i=1;i<=indbloc;i++) printf ("\n sed : %f",tra->vs[1][i]); */
  if ((resu=(struct resu*)malloc (sizeof(struct resu)))==NULL)
     {printf ("\n Error in allocation struct resu");
      exit(1);}
  moyo(tra,indbloc,npts,resu);
  il=strlen(nomsort);
  if (il<=6)
    {printf ("\n Le nom du fichier dat est mal interprete!!!!");
     exit(1);}
  il=il-6;
  
  for(lr=1;lr<3;lr++)
   {
     if (lr==1)    
        {
         strcpy(nombid,"mod.");
         strncat(nombid,nomsort,il);
         strcat(nombid,"HZ.dat.SPR");
         strcpy(nomin,"/home/alessia/pacific/scripts/mod.sprem.start.ray02");
         woutave(resu,nombid,nomin);
         fprintf (lssortie, "\n set Ficray%d = %s",ipath,nombid);
         strcpy(nombid,"mod.sour.");
         strncat(nombid,nomsort,il);
         strcat(nombid,"HZ.dat");
         fprintf (lssortie, "\n set Ficraysour%d = %s",ipath,nombid);
         strcpy(nomin2,"/home/alessia/pacific/scripts/mod.prem.iso.liss.ray02");
         woutsou(sour,nombid,nomin2);
        }
     else
        {
         strcpy(nombid,"mod.");
         strncat(nombid,nomsort,il);
         strcat(nombid,"HT.dat.SPR");
         strcpy(nomin,"/home/alessia/pacific/scripts/mod.sprem.start.lov02");
         woutave(resu,nombid,nomin);
         fprintf (lssortie, "\n set Ficlov%d = %s",ipath,nombid);
         strcpy(nombid,"mod.sour.");
         strncat(nombid,nomsort,il);
         strcat(nombid,"HT.dat");
         fprintf (lssortie, "\n set Ficlovsour%d = %s",ipath,nombid);
         strcpy(nomin2,"/home/alessia/pacific/scripts/mod.prem.iso.liss.lov02");
         woutsou(sour,nombid,nomin2);
        }
   }
  free(tra);
  free(resu);
  }
fclose (lssortie);
}   
/*#####write in the output files#######*/
/* write the average crustal model for the path in the output file.
   The output file is in the input format of saito94
   In this version the crustal structure is extended at depth with the upper mantle
   model PREM SMOOTHED with Q taken constant (Q=200).
   This upper mantle model (which must have 48 couches with its own crustal
   structure) and the input parameter for the saito97 pgm 
   are readen in the files mod.sprem.start.ray and mod.sprem.start.lov */ 

void woutave(resu,nomsort,nomin)
struct resu *resu;
char nomsort[100];
char nomin[100];
{
FILE *sortie;
FILE *entree;
char ligne[100];
int i,j,jtest;
int ncouch,nmod,nper;
float prof,profi,profint,profint2,epi,roi,vpi,vsi,xi,phi,eta,ai,ep;
float rho,lvp,lvs,lQ;
   
    
  if(!(sortie=fopen(nomsort,"w")))
     {printf ("%s\n : pb opening file",nomsort);
      exit(1);}
  if(!(entree=fopen(nomin,"r")))
       {printf ("\n unknown file : %s \n",nomin);
        exit(1);}
  
  printf ("\n reading file : %s",nomin);
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  sscanf(ligne,"%d ",&ncouch);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  /* the  average crustal structure is now written in the output file*/ 

  /*cas ou l'epaisseur de sediments est nulle : on ajoute 2 couches de
    croute sup de maniere a garder un nbre de couche constant*/
 
  if ( ((resu->ep[5]-0.021)<-1e-07) && (resu->ep[7]>=0.041))
       {
       resu->ep[5]=0.021;              /* "sediment" layer becomes part of upper crust*/
       resu->rho[4]=resu->rho[5]=resu->rho[6];
       resu->vp[4]=resu->vp[5]=resu->vp[6];
       resu->vs[4]=resu->vs[5]=resu->vs[6];
       resu->Q[4]=resu->Q[5]=resu->Q[6];
       printf("\n 0 resu->ep[7]=%f ",resu->ep[7]);
       resu->ep[7]=resu->ep[7]-0.021;
       printf("\n 1 resu->ep[7]=%f ",resu->ep[7]);
      
       }
   else if ( ((resu->ep[5]-0.021)<-1e-07) && (resu->ep[7]<0.041) )
       {
         printf ("\n No sediment and upper crust thickness < 40 m");
         exit(1);
       }

  /*cas ou l'epaisseur de glace est nulle : on ajoute 2 couches de
    croute sup de maniere a garder un nbre de couche constant*/
 
       if ( ((resu->ep[3]-0.021)<-1e-07) && (resu->ep[7]>=0.041) )
       {
       resu->ep[3]=resu->ep[5];     /* "ice" layer becomes "sediment" layer */
       resu->rho[2]=resu->rho[4];
       resu->rho[3]=resu->rho[5];
       resu->vp[2]=resu->vp[4];
       resu->vp[3]=resu->vp[5];
       resu->vs[2]=resu->vs[4];
       resu->vs[3]=resu->vs[5];
       resu->Q[2]=resu->Q[4];
       resu->Q[3]=resu->Q[5];
 
       resu->ep[5]=0.021;              /* "sediment" layer becomes part of upper crust*/       
       resu->rho[4]=resu->rho[5]=resu->rho[6];
       resu->vp[4]=resu->vp[5]=resu->vp[6];
       resu->vs[4]=resu->vs[5]=resu->vs[6];
       resu->Q[4]=resu->Q[5]=resu->Q[6];
 
       resu->ep[7]=resu->ep[7]-0.021;
       printf("\n 2 resu->ep[7]=%f ",resu->ep[7]);
       }
       else if ( ((resu->ep[3]-0.021)<-1e-07) && (resu->ep[7]<0.041) )
       {
         printf ("\n No Ice and upper crust thickness < 40 m");
         exit(1);
       }

  /*cas ou l'epaisseur d'eau est nulle : on ajoute 2 couches de
    croute sup de maniere a garder un nbre de couche constant*/
 
   if ( ((resu->ep[1]-0.021)<-1e-07) && (resu->ep[7]>=0.041))
       {
       resu->ep[1]=resu->ep[3];            /* "water" layer becomes "ice" layer */
       resu->rho[0]=resu->rho[2];
       resu->rho[1]=resu->rho[3];
       resu->vp[0]=resu->vp[2];
       resu->vp[1]=resu->vp[3];
       resu->vs[0]=resu->vs[2];
       resu->vs[1]=resu->vs[3];
       resu->Q[0]=resu->Q[2];
       resu->Q[1]=resu->Q[3];
 
       resu->ep[3]=resu->ep[5];            /* "ice" layer becomes sediment layer */
       resu->rho[2]=resu->rho[4];
       resu->rho[3]=resu->rho[5];
       resu->vp[2]=resu->vp[4];
       resu->vp[3]=resu->vp[5];
       resu->vs[2]=resu->vs[4];
       resu->vs[3]=resu->vs[5];
       resu->Q[2]=resu->Q[4];
       resu->Q[3]=resu->Q[5];
 
       resu->ep[5]=0.021;              /* "sediment" layer becomes part of upper crust*/      
       resu->rho[4]=resu->rho[5]=resu->rho[6];
       resu->vp[4]=resu->vp[5]=resu->vp[6];
       resu->vs[4]=resu->vs[5]=resu->vs[6];
       resu->Q[4]=resu->Q[5]=resu->Q[6];
 
       resu->ep[7]=resu->ep[7]-0.021;
       printf("\n 3 resu->ep[7]=%f ",resu->ep[7]);
       }
  else if ( ((resu->ep[1]-0.021)<-1e-07) && (resu->ep[7]<0.041) )
       {
         printf ("\n No Water and upper crust thickness < 40 m");
         exit(1);
       }

  prof=0.;
  for (j=0;j<9;j+=2)
    {  
     if(resu->Q[j]!=0) resu->Q[j]=1./600.;      
     fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f"
            ,prof,0.,resu->rho[j],resu->vp[j],resu->vs[j],0.,0.,0.);
     fprintf(sortie, "%7.4f\n",resu->Q[j]);
           
     ep=resu->ep[j+1]/2.;
     ep=(int)(100*ep)/100.;
     printf("\n\n resu->ep[%d]=%f ep=%f",j+1,resu->ep[j+1],ep);
     rho=( resu->rho[j] + resu->rho[j+1] ) / 2.;
     lvp=( resu->vp[j]  + resu->vp[j+1] ) / 2.;
     lvs=( resu->vs[j]  + resu->vs[j+1] ) / 2.;
 
     prof=prof+ep;
     fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f"
            ,prof,ep,rho,lvp,lvs,0.,0.,0.);
     fprintf(sortie, "%7.4f\n",resu->Q[j]);

     prof=prof+ep;
     if(resu->Q[j+1]!=0) resu->Q[j+1]=1./600.;      
     fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f"
            ,prof,ep,resu->rho[j+1],resu->vp[j+1],resu->vs[j+1],0.,0.,0.);
     fprintf(sortie, "%7.4f\n",resu->Q[j+1]);
    };
  
  /* the  average crustal structure is completed at depth with the upper mantle  model */
   
  jtest=0;
  for (j=1;j<=ncouch;j++)
     {
     fgets(ligne,100,entree);
     sscanf(ligne,"%f %f %f %f %f %f %f %f %f ",&profi,&epi,&roi,&vpi,&vsi,&xi,&phi,&eta,&ai);
     if((profi==40.) && (profi>prof))
         {
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                 ,prof,0.,roi,vpi,vsi,xi,phi,eta,0.005);
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                 ,profi,(profi-prof),roi,vpi,vsi,xi,phi,eta,0.005);
         jtest++;
         }
     else if((profi==50.) && (profi>(prof+0.02)) && (jtest==0))
         {
         profint=(profi+prof)/2.;
         profint=(int)(100*profint)/100.;  /*avoid problems when printing thickness 
                                             (prof+0.02) above for same reason      */
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,prof,0.,roi,vpi,vsi,xi,phi,eta,0.005);
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,profint,(profint-prof),roi,vpi,vsi,xi,phi,eta,0.005);
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,profi,(profi-profint),roi,vpi,vsi,xi,phi,eta,0.005);
         jtest++;
         }
     else if((profi==75.) && (profi>(prof+0.04)) && (jtest==0))
         {
         profint=(profi+prof)/2.;
         profint=(int)(100*profint)/100.; /*avoid problems when printing thickness 
                                            (prof+0.02) above for same reason      */
         profint2=(profint+prof)/2.;
         profint2=(int)(100*profint2)/100.;
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,prof,0.,roi,vpi,vsi,xi,phi,eta,0.005);
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,profint2,(profint2-prof),roi,vpi,vsi,xi,phi,eta,0.005);
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,profint ,(profint-profint2),roi,vpi,vsi,xi,phi,eta,0.005);
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,profi,(profi-profint),roi,vpi,vsi,xi,phi,eta,0.005);
         jtest++;
         }
     else if((profi>=75.) && (jtest==0))
         {
         printf ("\n average crustal thickness of path %s >=75 km? \nMay be an error....", nomin);
         exit(1);
         }
     else if (!(jtest==0.))
         {
         fprintf(sortie, "%7.2f %7.2f %6.3f %6.3f %6.3f%7.4f%7.4f%7.4f%7.4f\n"
                ,profi,epi,roi,vpi,vsi,xi,phi,eta,0.005); 
         }
     }
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  sscanf(ligne,"%d ",&nmod);
  fprintf(sortie,"   1                   nmodes\n");
/*fputs(ligne,sortie);*/
  
  for (i=1;i<=nmod;i++)
    {
     fgets(ligne,100,entree);
     sscanf(ligne,"%d ",&nper);
     fputs(ligne,sortie);
     for (j=1;j<=nper;j++)
       {
       fgets(ligne,100,entree);
       fputs(ligne,sortie);
       }
    }
  fclose(entree);
  fclose(sortie);
}
/*#####write in the output files#######*/
/* write the source mantel model for the path in the output file.
   After 350 km the source model is extended at depth with a smoothed
   prem model.
   The output file is in the input format of saito94
   structure and the input parameter for the saito94 pgm
   are readen in the files mod.prem.iso.liss.ray  and 
   mod.prem.iso.liss.lov.    */
 
void woutsou(sour,nomsort,nomin)
struct source *sour;
char nomsort[100];
char nomin[100];
{
FILE *sortie;
FILE *entree;
char ligne[100];
int i,j,jtest,istart;
int ncouch,nmod,nper;
float prof,profi,profint,profint2,profint3,epi,roi,vpi,vsi,bid,ai,ep;
float rho,lvp,lvs,lQ;
  
  if(!(sortie=fopen(nomsort,"w")))
     {printf ("%s\n : pb opening file",nomsort);
      exit(1);}  
  if(!(entree=fopen(nomin,"r")))
       {printf ("\n unknown file : %s \n",nomin);
        exit(1);}
 
  printf ("\n reading file : %s",nomin);
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
 
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
 
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
 
  fgets(ligne,100,entree);
  sscanf(ligne,"%d ",&ncouch);
  fputs(ligne,sortie);
 
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
 
  fgets(ligne,100,entree);
  fputs(ligne,sortie);

  /* the  average crustal structure is now written in the output file*/

  /*cas ou l'epaisseur de sediments est nulle : on ajoute 2 couches de
    croute sup de maniere a garder un nbre de couche constant*/

  if ( ((sour->ep[5]-0.021)<-1e-07) && (sour->ep[7]>=0.041))
       {
       sour->ep[5]=0.021;              /* "sediment" layer becomes part of upper crust*/
       sour->rho[4]=sour->rho[5]=sour->rho[6];
       sour->vp[4]=sour->vp[5]=sour->vp[6];
       sour->vs[4]=sour->vs[5]=sour->vs[6];
       sour->Q[4]=sour->Q[5]=sour->Q[6]; 
       sour->ep[7]=sour->ep[7]-0.021;
       }
   else if ( ((sour->ep[5]-0.021)<-1e-07) && (sour->ep[7]<0.041) )
       {
         printf ("\n No sediment and upper crust thickness < 40 m");
         exit(1);
       }
  /*cas ou l'epaisseur de glace est nulle : on ajoute 2 couches de
    croute sup de maniere a garder un nbre de couche constant*/
 
       if ( ((sour->ep[3]-0.021)<-1e-07) && (sour->ep[7]>=0.041))
       {
       sour->ep[3]= sour->ep[5];            /* "ice" layer becomes "sediment" layer */
       sour->rho[2]=sour->rho[4];
       sour->rho[3]=sour->rho[5];
       sour->vp[2]=sour->vp[4];
       sour->vp[3]=sour->vp[5];
       sour->vs[2]=sour->vs[4];
       sour->vs[3]=sour->vs[5];
       sour->Q[2]=sour->Q[4];
       sour->Q[3]=sour->Q[5];
 
       sour->ep[5]=0.021;              /* "sediment" layer becomes part of upper crust*/
       sour->rho[4]=sour->rho[5]=sour->rho[6];
       sour->vp[4]=sour->vp[5]=sour->vp[6];
       sour->vs[4]=sour->vs[5]=sour->vs[6];
       sour->Q[4]=sour->Q[5]=sour->Q[6];
 
       sour->ep[7]=sour->ep[7]-0.021;
       }
  else if ( ((sour->ep[3]-0.021)<-1e-07) && (sour->ep[7]<0.041) )
       {
         printf ("\n No Ice and upper crust thickness < 40 m");
         exit(1);
       }

  /*cas ou l'epaisseur d'eau est nulle : on ajoute 2 couches de
    croute sup de maniere a garder un nbre de couche constant*/

   if ( ((sour->ep[1]-0.021)<-1e-07) && (sour->ep[7]>=0.041))
       {
       sour->ep[1]= sour->ep[3];            /* "water" layer becomes "ice" layer */
       sour->rho[0]=sour->rho[2];
       sour->rho[1]=sour->rho[3];
       sour->vp[0]=sour->vp[2];
       sour->vp[1]=sour->vp[3];
       sour->vs[0]=sour->vs[2];
       sour->vs[1]=sour->vs[3];
       sour->Q[0]=sour->Q[2];
       sour->Q[1]=sour->Q[3];

       sour->ep[3]= sour->ep[5];            /* "ice" layer becomes sediment layer */
       sour->rho[2]=sour->rho[4];
       sour->rho[3]=sour->rho[5];
       sour->vp[2]=sour->vp[4];
       sour->vp[3]=sour->vp[5];
       sour->vs[2]=sour->vs[4];
       sour->vs[3]=sour->vs[5];
       sour->Q[2]=sour->Q[4];
       sour->Q[3]=sour->Q[5];
 
       sour->ep[5]=0.021;              /* "sediment" layer becomes part of upper crust*/
       sour->rho[4]=sour->rho[5]=sour->rho[6];
       sour->vp[4]=sour->vp[5]=sour->vp[6];
       sour->vs[4]=sour->vs[5]=sour->vs[6];
       sour->Q[4]=sour->Q[5]=sour->Q[6];
 
       sour->ep[7]=sour->ep[7]-0.021;
       }
  else if ( ((sour->ep[1]-0.021)<-1e-07) && (sour->ep[7]<0.041) )
       {
         printf ("\n No Water and upper crust thickness < 40 m");
         exit(1);
       }

  prof=0; 
  for (j=0;j<9;j+=2)
    {
     if(sour->vs[j]==0)  lvs=0;
     else                lvs=1/sour->vs[j];
     if(sour->vp[j]==0)  lvp=0;
     else                lvp=1/sour->vp[j];
     fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f"
            ,prof,0.,sour->rho[j],lvp,lvs,0.,0.,0.);
     fprintf(sortie, " %6.4f\n",sour->Q[j]);

     ep=sour->ep[j+1]/2.;
     ep=(int)(100*ep)/100.;
     rho=( sour->rho[j] + sour->rho[j+1] ) / 2.;
     if( !(sour->vp[j]==0)  && !(sour->vp[j+1]==0) )
           lvp=( (1/sour->vp[j]) + (1/sour->vp[j+1]) ) / 2.;
     else  lvp=0.;
     if( !(sour->vs[j]==0)  && !(sour->vs[j+1]==0) )
           lvs=( (1/sour->vs[j]) + (1/sour->vs[j+1]) ) / 2.;
     else  lvs=0.;
     lQ=( sour->Q[j] + sour->Q[j+1] ) / 2.;
     prof=prof+ep;
     fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f"
            ,prof,ep,rho,lvp,lvs,0.,0.,0.);
     fprintf(sortie, " %6.4f\n",lQ);

     if(sour->vs[j+1]==0)  lvs=0;
     else                lvs=1/sour->vs[j+1];
     if(sour->vp[j+1]==0)  lvp=0;
     else                lvp=1/sour->vp[j+1];
     prof=prof+ep;
     fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f"
            ,prof,ep,sour->rho[j+1],lvp,lvs,0.,0.,0.);
     fprintf(sortie, " %6.4f\n",sour->Q[j+1]);
    };

  /* the  source crustal structure is completed at depth with the upper mantle  crustal
     structure of 3SMAC until a depth of 350 km. Below 350 km we use a smoothed prem model */
  /* we add the upper mantle structure of 3SMAC for the first 350 km */

     if((50.-prof)>0.)
         {
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,prof,0.,sour->rho[10],1/sour->vp[10],1/sour->vs[10],0.0,0.0,0.0,sour->Q[10]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,50.,50.-prof,sour->rho[10],1/sour->vp[10],1/sour->vs[10],0.0,0.0,0.0,sour->Q[10]);
         istart=11;
         profint=50;
         }
     else if ((60-(prof+0.02))>0)
         {
         profint =(60+prof)/2;
         profint=(int)(100*profint)/100.;  /* avoid problems when printing thickness
                                              same for (prof+0.02) above */ 
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,prof,0.,sour->rho[11],1/sour->vp[11],1/sour->vs[11],0.0,0.0,0.0,sour->Q[11]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profint,profint-prof,sour->rho[11],1/sour->vp[11],1/sour->vs[11],0.0,0.0,0.0,sour->Q[11]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,60.,60.-profint,sour->rho[11],1/sour->vp[11],1/sour->vs[11],0.0,0.0,0.0,sour->Q[11]);
         istart=12;
         profint=60;
         }
     else if ((70-(prof+0.04))>0)
         {
         profint =(70+prof)/2;
         profint=(int)(100*profint)/100.;  /* avoid problems when printing thickness
                                              same for (prof+0.02) above */ 
         profint2 =(70+profint)/2;
         profint2=(int)(100*profint2)/100.; 
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,prof,0.,sour->rho[12],1/sour->vp[12],1/sour->vs[12],0.0,0.0,0.0,sour->Q[12]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profint,profint-prof,sour->rho[12],1/sour->vp[12],1/sour->vs[12],0.0,0.0,0.0,sour->Q[12]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profint2,profint2-profint,sour->rho[12],1/sour->vp[12],1/sour->vs[12],0.0,0.0,0.0,sour->Q[12]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,70.,70-profint2,sour->rho[12],1/sour->vp[12],1/sour->vs[12],0.0,0.0,0.0,sour->Q[12]);
         istart=13;
         profint=70;
         }
     else if ((80-(prof+0.08))>0)
        {
         profint =(80+prof)/2;           /* avoid problems when printing thickness
                                            same for (prof+0.08) above */
         profint=(int)(100*profint)/100.; 
         profint2 =(80+profint)/2;
         profint2=(int)(100*profint2)/100.; 
         profint3 =(80+profint2)/2;
         profint3=(int)(100*profint3)/100.; 
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,prof,0.,sour->rho[13],1/sour->vp[13],1/sour->vs[13],0.0,0.0,0.0,sour->Q[13]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profint,profint-prof,sour->rho[13],1/sour->vp[13],1/sour->vs[13],0.0,0.0,0.0,sour->Q[13]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profint2,profint2-profint,sour->rho[13],1/sour->vp[13],1/sour->vs[13],0.0,0.0,0.0,sour->Q[13]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profint3,profint3-profint,sour->rho[13],1/sour->vp[13],1/sour->vs[13],0.0,0.0,0.0,sour->Q[13]);
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,80.,80-profint3,sour->rho[13],1/sour->vp[13],1/sour->vs[13],0.0,0.0,0.0,sour->Q[13]);
         istart=14;
         profint=80;
         }
  for (j=istart;j<24;j++)
         {
         if      (j==10) prof=50.;
         else if (j==11) prof=60.;
         else if (j==12) prof=70.;
         else if (j==13) prof=80.;
         else if (j==14) prof=90.;
         else if (j==15) prof=100.;
         else if (j==16) prof=115.;
         else if (j==17) prof=130.;
         else if (j==18) prof=150.;
         else if (j==19) prof=175.;
         else if (j==20) prof=200.;
         else if (j==21) prof=250.;
         else if (j==22) prof=300.;
         else if (j==23) prof=350.;
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,prof,prof-profint,sour->rho[j],1/sour->vp[j],1/sour->vs[j],0.0,0.0,0.0,sour->Q[j]);
         profint=prof;
         }
/* we add the smoothed prem model */
  jtest=0;
  for (j=1;j<=ncouch;j++)
     {
     fgets(ligne,100,entree);
     sscanf(ligne,"%f %f %f %f %f %f %f %f %f ",&profi,&epi,&roi,&vpi,&vsi,&bid,&bid,&bid,&ai); 
     
     if((profi>350) && (jtest==0)) 
        {
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profi,profi-350.,roi,vpi,vsi,bid,bid,bid,ai); 
         jtest=1;
        }
     else if((profi>=350) && (jtest==1))
        { 
         fprintf(sortie, "%7.2f %6.2f %6.3f %6.3f %6.3f %6.4f %6.4f %6.4f %6.4f\n"
                ,profi,epi,roi,vpi,vsi,bid,bid,bid,ai);
        }
     }    
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  fputs(ligne,sortie);
  
  fgets(ligne,100,entree);
  sscanf(ligne,"%d ",&nmod);
  fprintf(sortie,"   1                   nmodes\n");
/*fputs(ligne,sortie);*/
  
  for (i=1;i<=nmod;i++)
    {
     fgets(ligne,100,entree);
     sscanf(ligne,"%d ",&nper);
     fputs(ligne,sortie);
     for (j=1;j<=nper;j++)
       { 
       fgets(ligne,100,entree);
       fputs(ligne,sortie);
       } 
    }
  fclose(entree);
  fclose(sortie);
}
/*#####average calculation#######*/
void moyo(toto,nbloc,npts,resu)
int npts[400];
struct layer *toto;
struct resu *resu;
int nbloc;
{
int i,j,some;
int somvp[10],somvs[10],somrho[10],somQ[10];
float la,lo;
 
 resu->ep[0]= 0.;                                   
 resu->ep[1]= 0.;                                   
 resu->ep[2]= 0.;                                   
 resu->ep[3]= 0.;                                   
 resu->ep[4]= 0.;                                   
 resu->ep[5]= 0.;                                   
 resu->ep[6]= 0.;                                   
 resu->ep[7]= 0.;                                   
 resu->ep[8]= 0.;                                   
 resu->ep[9]= 0.;                                   
 some=0.;                                   
 for (i=1; i<=nbloc; i++)
      {
        resu->ep[1]= resu->ep[1] + npts[i]*toto->epe[i];
        resu->ep[3]= resu->ep[3] + npts[i]*toto->epg[i];
        resu->ep[5]= resu->ep[5] + npts[i]*toto->eps[i];
        resu->ep[7]= resu->ep[7] + npts[i]*toto->epcs[i];
        resu->ep[9]= resu->ep[9] + npts[i]*toto->epci[i];
        some=some+npts[i];
      }
    resu->ep[1]=resu->ep[1]/(float)some;
    resu->ep[3]=resu->ep[3]/(float)some;
    resu->ep[5]=resu->ep[5]/(float)some;
    resu->ep[7]=resu->ep[7]/(float)some;
    resu->ep[9]=resu->ep[9]/(float)some;
 
 /* avoid problems when printing ep*/ 
    resu->ep[1]=(int)(100*resu->ep[1])/100.;
    resu->ep[3]=(int)(100*resu->ep[3])/100.;
    resu->ep[5]=(int)(100*resu->ep[5])/100.;
    resu->ep[7]=(int)(100*resu->ep[7])/100.;
    resu->ep[9]=(int)(100*resu->ep[9])/100.;
   

for (j=0; j<10; j++)
{
 somvp[j]=0;
 somvs[j]=0;
 somrho[j]=0;
 somQ[j]=0;
 resu->vp[j]=0.;                                    
 resu->vs[j]=0.;
 resu->rho[j]=0.;
 resu->Q[j]=0.;
     /*   printf ("\n NBLOC=%d \n",nbloc); */
    for (i=1; i<=nbloc; i++)
      { 
        resu->vp[j]=resu->vp[j] + npts[i]*toto->vp[j][i];
/*      printf ("\n M : resu->vp[j]=%f npts[%d]=%d toto->vp[j][i]=%f ",resu->vp[j],i,npts[i],toto->vp[j][i]); */
        if (!(toto->vp[j][i]==0.)) somvp[j]=somvp[j] + npts[i];

        resu->vs[j]=resu->vs[j] + npts[i]*toto->vs[j][i];
        if (!(toto->vs[j][i]==0.)) somvs[j]=somvs[j] + npts[i];

        resu->rho[j]=resu->rho[j] + npts[i]*toto->rho[j][i];
        if (!(toto->rho[j][i]==0.)) somrho[j]=somrho[j] + npts[i];

        resu->Q[j]=resu->Q[j] + npts[i]*toto->Q[j][i];
        if (!(toto->Q[j][i]==0.)) somQ[j]=somQ[j] + npts[i];
      }

    if (!(resu->vp[j]==0.)) resu->vp[j]=1/(resu->vp[j]/(float)somvp[j]);
    else if (j>5)                                                  
    {printf ("\n MOYO : resu->vp[j]==%f at layer %d \n", resu->vp[j],j);
     exit(1);}

    if (!(resu->vs[j]==0.)) resu->vs[j]=1/(resu->vs[j]/(float)somvs[j]);
    else if (j>5)                                                  
    {printf ("\n MOYO : resu->vs[j]==0. at layer %d \n",j);
     exit(1);}

    if (!(resu->Q[j]==0.)) resu->Q[j]=resu->Q[j]/(float)somQ[j];
    else if (j>5)                                                   
    {printf ("\n MOYO : resu->Q[j]==0. at layer %d \n",j);
     exit(1);}

    if (!(resu->rho[j]==0.)) resu->rho[j]=resu->rho[j]/(float)somrho[j];
    else if (j>5)                                                  
    {printf ("\n MOYO : resu->rho[j]==0. at layer %d \n",j);
     exit(1);}
 }
}
/*#####Affectation of the 3SMAC structure###########*/
void affstructsmac(smac,mant)
struct smac *smac;
struct mant *mant;
{
int i,j,jj,kk,bb,comptt,cvp,cvs,crho,cq;
char nomfich[100];
char nomfich2[100];
float incrst[180][90];
int indep[180][90];
float **ti;
float **ta;
float **tb;
float **tc;
float **td;
float **tcr;


/*  Read the ice data */
  sprintf (nomfich,"brut/dataice");
  ti=lecsmac(nomfich); 
for(i=0;i<180;i++)
     {
     for(j=0;j<90;j++)
          {
           smac->ice[i][j]=ti[i][j];
           indep[i][j]=0;
      /*   printf("\n  smac->ice[%d][%d] = %f ",i,j,smac->ice[i][j]); */
          }
     }
free(ti);


/*layer thickness  */

for(jj=1;jj<7;jj+=2) /*jj+=2 equivaut jj=jj+2 */
{
  sprintf (nomfich,"%s.%d","radi/radi",jj);
  ta=lecsmac(nomfich); 
  sprintf (nomfich,"%s.%d","radi/radi",jj+1);
  tb=lecsmac(nomfich);
  for(i=0;i<180;i++)
     {
     for(j=0;j<90;j++)
           {
            if( (smac->ice[i][j]>0) && (jj==1) )
               {
                smac->ep[indep[i][j]][i][j]=0;
                indep[i][j]++;
                smac->ep[indep[i][j]][i][j]=ta[i][j]-tb[i][j];
                indep[i][j]++;
               }
            else if ( (smac->ice[i][j]==0) && (jj==1) )
               {
                smac->ep[indep[i][j]][i][j]=ta[i][j]-tb[i][j]; 
                indep[i][j]++;  
                smac->ep[indep[i][j]][i][j]=0;
                indep[i][j]++;
               }
            else if (jj>1)
               {
                smac->ep[indep[i][j]][i][j]=ta[i][j]-tb[i][j]; 
                indep[i][j]++;
               }

            ta[i][j]=0.;
            tb[i][j]=0.;
            if (indep[i][j]>4)
               {printf("\n indep[i][j] = %d ",indep[i][j]);
                exit(1);}
           }
     }
free(ta);
free(tb);
} 

/*lower crust thickness :  the file datacroute
  contain only the crustal thickness and does not include sediments (ok)
  To check this look in the pgm make_radi.f given by Nataf et Ricard
  with their model.
  simple method */

sprintf (nomfich,"%s","brut/datacroute");
tcr=lecsmac(nomfich);
  for(i=0;i<180;i++)
     {   
     for(j=0;j<90;j++)
           {
            smac->ep[4][i][j]=(tcr[i][j]/1000.)-smac->ep[3][i][j];
            incrst[i][j]=0;
           }
     }   
free(tcr);

/*lecture des fichiers crust.X qui permettront de determiner
la profondeur ou lire les parametres elastiques (affectation
de incrst[i][j] */

for(jj=9;jj<=15;jj++)
  {
 sprintf (nomfich2,"%s.%d","radi/crust",jj);
 tcr=lecsmac(nomfich2);
  for(i=0;i<180;i++)  
     {   
     for(j=0;j<90;j++)   
       {  
       if((tcr[i][j]==0.)&&(incrst[i][j]==0)) incrst[i][j]=jj-1;
       } 
     }
  }
free(tcr);



/*elastic parameters */
/*  jj=1 and 2 are for the water or ice layer*/

comptt=0;bb=0;
for(jj=1;jj<27;jj++)
 {
  cvp=cvs=crho=cq=0;
  sprintf (nomfich,"%s.%d","para/VP",jj);
  printf ("\n %s.%d","para/VP",jj);
  ta=lecsmac(nomfich);    
  sprintf (nomfich,"%s.%d","para/VS",jj);
  tb=lecsmac(nomfich);    
  sprintf (nomfich,"%s.%d","para/RHO",jj);
  tc=lecsmac(nomfich);    
  sprintf (nomfich,"%s.%d","para/LnQMU",jj);
  td=lecsmac(nomfich);   
  for(i=0;i<180;i++) 
     {   
     for(j=0;j<90;j++)
       {
       if( (smac->ice[i][j]>0) && (jj<=2) )
        {
         smac->vp[jj-1][i][j]=1./1.450;                            /*  water layer */
         smac->vs[jj-1][i][j]=0.;
         smac->rho[jj-1][i][j]=1.020;
         smac->Qmu[jj-1][i][j]=0.0017;
         kk=(jj-1)+2;

         if(ta[i][j]>=99999) {smac->vp[kk][i][j]=0.;  cvp++;}      /*  ice layer  */
         else                smac->vp[kk][i][j]=1/ta[i][j];
 
         if(tb[i][j]>=99999) {smac->vs[kk][i][j]=0.; cvs++;} 
         else                smac->vs[kk][i][j]=1/tb[i][j];
 
         if(tc[i][j]>=99999) {smac->rho[kk][i][j]=0.; crho++;}
         else                smac->rho[kk][i][j]=tc[i][j];
 
         if(td[i][j]>=99999) {smac->Qmu[kk][i][j]=0.; cq++;}
         else                smac->Qmu[kk][i][j]=1./exp(td[i][j]);
        } 
       else if ( (smac->ice[i][j]==0) && (jj<=2) )
        {
         if(ta[i][j]>=99999) {smac->vp[jj-1][i][j]=0.; cvp++;}       /* water layer */
         else                smac->vp[jj-1][i][j]=1/ta[i][j];
 
         if(tb[i][j]>=99999) {smac->vs[jj-1][i][j]=0.; cvs++;}
         else                smac->vs[jj-1][i][j]=0.;
 
         if(tc[i][j]>=99999) {smac->rho[jj-1][i][j]=0.; crho++;}
         else                smac->rho[jj-1][i][j]=tc[i][j];
 
         if(td[i][j]>=99999) {smac->Qmu[jj-1][i][j]=0.; cq++;}
         else                smac->Qmu[jj-1][i][j]=0.;  /* 1/Qmu set to zero because
                                                           Qmu does not exist in water
                                                           Qmu =0 is preferred to infinity to
                                                           avoid numerical problems in saito02
                                                           Any value of Qmu in water should anyway
                                                           not affect the result                   */
         kk=(jj-1)+2;
         smac->vp[kk][i][j]=0.;                           /*  ice layer */
         smac->vs[kk][i][j]=0.;
         smac->rho[kk][i][j]=0;
         smac->Qmu[kk][i][j]=0;
        }
       else if ((jj>2)&&(jj<8))
        {
         if(ta[i][j]>=99999) {smac->vp[jj+1][i][j]=0.; cvp++;}
         else                smac->vp[jj+1][i][j]=1/ta[i][j];

         if(tb[i][j]>=99999) {smac->vs[jj+1][i][j]=0.; cvs++;}
         else                smac->vs[jj+1][i][j]=1/tb[i][j];
         
         if(tc[i][j]>=99999) {smac->rho[jj+1][i][j]=0.; crho++;}
         else                smac->rho[jj+1][i][j]=tc[i][j];
         
         if(td[i][j]>=99999) {smac->Qmu[jj+1][i][j]=0.; cq++;}
         else                smac->Qmu[jj+1][i][j]=1./exp(td[i][j]);
        }
       else if ((jj>=8)&&(jj<16)&&(jj==incrst[i][j]))
        {
         if(ta[i][j]>=99999) {smac->vp[9][i][j]=0.; bb=1;}
         else                smac->vp[9][i][j]=1/ta[i][j];

         if(tb[i][j]>=99999) {smac->vs[9][i][j]=0.; bb=1;}
         else                smac->vs[9][i][j]=1/tb[i][j];

         if(tc[i][j]>=99999) {smac->rho[9][i][j]=0.; bb=1;}
         else                smac->rho[9][i][j]=tc[i][j];

         if(td[i][j]>=99999) {smac->Qmu[9][i][j]=0.; bb=1;}
         else                smac->Qmu[9][i][j]=1./exp(td[i][j]);
         if(bb==0)comptt++;
        }
       if(jj>11)
        {
         if(ta[i][j]>=99999) mant->vp[jj-12][i][j]=0.;
         else                mant->vp[jj-12][i][j]=1/ta[i][j];

         if(tb[i][j]>=99999) mant->vs[jj-12][i][j]=0.;
         else                mant->vs[jj-12][i][j]=1/tb[i][j];

         if(tc[i][j]>=99999) mant->rho[jj-12][i][j]=0.;
         else                mant->rho[jj-12][i][j]=tc[i][j];

         if(td[i][j]>=99999) mant->Qmu[jj-12][i][j]=0.;
         else                mant->Qmu[jj-12][i][j]=1./exp(td[i][j]);
        }
       }
     }   
  if (jj<8)
    {printf ("\n Number of blocks for which a parameter is not determined : ");
     printf ("\n Layer %d  VP : %d;  VS : %d;  RHO : %d;  Qmu : %d \n",jj,cvp,cvs,crho,cq);}
  free(ta);
  free(tb);
  free(tc);
  free(td);   
 }
  if (comptt-16200)
    {printf ("\n elastic parameters just above Moho are not everywhere determined?");
     printf ("\n comptt = %d",comptt);
     exit(1);} 
}
/*#####reading the 3SMAC model ##############*/
float ** lecsmac(nomfich)
char *nomfich;
{
int i, j;
char nom[200];
char bid[200];
char ligne1[300];
char ligne2[300];
char ligne3[300];
float **par;
FILE *entree;

par=allotab(180,90);
strcpy(nom,"/home/alessia/pacific/3-SMAC/");
strcat(nom,nomfich);
if(!(entree=fopen(nom,"r")))
     {printf ("\n unknown file : %s \n",nom);
      exit(1);}
   printf ("\n reading file : %s",nom);   
fgets(ligne1,300,entree);
fgets(ligne2,300,entree);
fgets(ligne3,300,entree);
for(i=0;i<180;i++)
     {
/*   if(i==179) 
            printf("\t lecsmac: i=179\n");  */
     for(j=0;j<90;j++)
          {
          fscanf (entree, "%e", &par[i][j]);
          }
     }
fclose(entree);
return(par);
}
/*#### memory allocation for matrice and vector #######*/
float** allotab(l,c)
int l,c;
{
float **mat;
int i;
if ((mat=(float**)malloc(l*sizeof(float*)))==NULL)
   {printf ("Error in matrix allocation 1");
    exit(1);}
    for (i=0;i<l;i++)
     {   
     if ((mat[i]=(float*)malloc(c*sizeof(float)))==NULL)
       {printf ("Error in matrix allocation 2");
        exit(1);}
     }
return(mat);   
}
/*######distance between 2 points on great circle######*/
/*The argument late,lone,lats,lons are supposed to be given 
  in radians. late and lats are colatitude (colatitude)*/

float  del(late,lone,lats,lons)
float late,lone,lats,lons;
{
  float resu,ca,cb,sa,sb;
 
 
  ca=cos(late);
  cb=cos(lats);
  sa=sin(late);
  sb=sin(lats);
  resu=ca*cb+(sa*sb)*cos(lone-lons);

  if(resu>1.)resu=1.;
  if(resu<-1.)resu=-1.;
  resu=acos(resu);
  return (resu);
}

/* ##########great cicle azimuth########################## */
/*The argument late,lone,lats,lons,delt are supposed to be given
  in radians. late and lats are latitude (not colatitude)*/

float azimuth(late,lone,lats,lons,delt)
float late,lone,lats,lons,delt;
{
#define Epsilon 0.00000001
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

/* #################################### */
/* computation of lati and loni for a point i located on the great circle 
   at a distance delt of the epicenter(late,lone). The argument late,lone,
   delt,azes are supposed to be in radians but the results lati and loni are
   converted in degrees.late is colatitude (not latitude)*/

void coord(late,lone,delt,azes)
float late,lone,delt,azes;
{
float ca,cb,cc,sa,sb,sc,cosA,C,cosC,val,lati1;

 /*calcul de la latitude */

  cosA=cos(azes);
  cb=cos(late);
  sb=sin(late);
  cc=cos(delt);
  sc=sin(delt);

  ca=cb*cc+sb*sc*cosA;
  if(ca>1.)ca=1.;
  if(ca<-1.)ca=-1.;
     lati1=(pi/2)-acos(ca);
     lati=(1/0.993277)*tan(lati1);
/*   lati=tan(lati1); */
     lati=180/pi*atan(lati);

 /*calcul de la longitude*/

   sa=sin((pi/2)-(lati1));

   cosC=(cc-ca*cb)/(sa*sb);
   if(cosC>1.)cosC=1.;
   if(cosC<-1.)cosC=-1.;
   C=acos(cosC);
   if(azes>=0 && azes<=pi) loni=C+lone;
   if(azes>pi && azes<2*pi) loni=lone-C;
   loni=loni*180/pi;
  }
