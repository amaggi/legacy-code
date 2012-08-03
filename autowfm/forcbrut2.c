#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
/*Theses values have to fit the one in nn.param
np_max is the maximum number of points
nd_max is the maximum number of dimensions (3)
nt_max is the maximum number of tetrahedra  */
#define np_max 65000
#define nd_max 3
#define nt_max 260000
#define nwork3d 650000 

/* Estimate the quality criterion for each voromoi cell  
In input of this pgm we need  2 files:
LISTNODES contains lon,lat for each nodes :
eg :
 143.178604    -13.959400   
 144.485199    -17.364401  
intomodesVs.GLOB contains  coord epicenter and station
*/

/* ##### ROUTINES ########### */
double ** allotabd();
int ** allotabi();
double  del();
double  azimuth();
void coord();
extern void nnspherericinit_(int *,int *,int *,double *,double *,int *);
extern void findnodesph_(int *,double *,int *,int *,double *,double *,int *);

/* ##### GLOBAL CONSTANTS */
double pi = 3.14159265358979323846264338;
double radtodeg; 
double degtorad;
double pis2; 

int main() {
  char comp,comp2;
  char ligne[400],input[100];
  char **qualitaz;
  char nomdatray[80],bids[13];
  char nombid[13],nombid2[13];
  char nomsta[800][5];
  int nt,jj,ii,ij,ijmem,ik,ilen,ichar;
  int intrala,ista,nsta,icompt;
  int ipath,npaths,nnodes,irun,nbid;
  int    *itest;
  int    *nopathpernode; 
  float ela,elo,sla,slo,lonod,lanod;

  float *stalat,*stalon;
  double *qualit;
  double *length,*azcell;
  double dst,dpas;
  double *coelat,*coslat;
  double *elon,*slon;
  double *colatnode,*lonnode;
  double distpath,azim,azimem;
  double colati,loni;
  double azes,delta;
  FILE *indat;
  FILE *inodes;
  FILE *output;
  /*#################### Parameters for NN routines######################  */
  int dimxd,iwalk,idebug,iwrite,mode_del,node;
  double xlon,xlat;
  double xd[2];

  /* set values of global constants */
  radtodeg = 180./pi;
  degtorad = pi/180.;
  pis2 = pi/2.;

  iwalk=0;
  dimxd=2;
  node = 2;

  idebug = 1;   /*idebug = 0 do not perform debug checks
                  idebug = 1 perform debug checks       */
  iwrite = 1;   /*iwrite = 1 write out vertices and neighbours
                  iwrite = 0 do not write information   */
  mode_del = 0; /*mode_del = 1 means read in delaunay
                  mode_del = 0 means calculate delaunay */

  /*INPUT FILES */
  printf("\n Name input file for forcbrut?  ");
  scanf("%s",input);
  if(!(inodes=fopen(input,"r"))) {
    printf ("\n forcbrut : pb opening input file %s",input);
    exit(1);
  }

  fgets(ligne,400,inodes);
  sscanf(ligne,"%d %d %d",&nnodes,&irun,&nbid);
  printf("\n%d nodes\n",nnodes);

  /* allocate memory for nodes*/
  if ((nopathpernode=(int*)malloc(nnodes*sizeof(int)))==NULL) {
    printf ("Error in vector allocation nopathpernode");
    exit(1);
  }
  if ((itest=(int*)malloc(nnodes*sizeof(int)))==NULL){
    printf ("Error in vector allocation itest");
    exit(1);
  }

  if ((lonnode=(double*)malloc(nnodes*sizeof(double)))==NULL)
    {printf ("Error in vector allocation lonnode");
     exit(1);}
  if ((colatnode=(double*)malloc(nnodes*sizeof(double)))==NULL)
    {printf ("Error in vector allocation colatnode");
     exit(1);}
  if ((length=(double*)malloc(nnodes*sizeof(double)))==NULL)
      {printf ("Error in vector allocation length");
        exit(1);}
  if ((azcell=(double*)malloc(nnodes*sizeof(double)))==NULL)
      {printf ("Error in vector allocation azcell");
        exit(1);}
  if ((qualit=(double*)malloc(nnodes*sizeof(double)))==NULL)
    {printf ("Error in vector allocation qualit");
     exit(1);}

  if ((qualitaz=(char**)malloc(nnodes*sizeof(char*)))==NULL) {
    printf ("Error 1  in matrix allocation qualitaz");
    exit(1);
  }
  for (ii=0;ii<nnodes;ii++) {
    if ((qualitaz[ii]=(char*)malloc(5*sizeof(char)))==NULL) {
      printf ("Error 2  in matrix allocation qualitaz");
      exit(1);
    }
  }
  /* end of memory allocation for nodes */


  /* read the nodes */
  for (ii=0;ii<nnodes;ii++) {
    fgets(ligne,400,inodes);
    sscanf(ligne,"%f %f",&lonod,&lanod);

    lonnode[ii]=(double)lonod;
    colatnode[ii]=(double)lanod;
    colatnode[ii]=pis2-(atan(0.993277*tan(colatnode[ii]*degtorad)));
    lonnode[ii]=lonnode[ii]*degtorad; 
  }
  fclose(inodes);

  nnspherericinit_(&idebug,&iwrite,&mode_del,colatnode,lonnode,&nnodes);
  printf(" End nnsphereric \n");
  /* end of nodes input */


  /*intomodes input */
  printf ("\n Lecture intomodesVs.GLOB"); 
  if(!(indat=fopen("intomodesVs.GLOB","r"))){
    printf ("\n pb opening file intomodesVs.GLOB \n");
    exit(1);
  }
  fgets(ligne,400,indat);
  sscanf(ligne,"%d %d ",&npaths);

  /* allocate memory for paths */
  if ((coelat=(double*)malloc(npaths*sizeof(double)))==NULL)
      {printf ("Error in vector allocation coelat");
        exit(1);}
  if ((elon=(double*)malloc(npaths*sizeof(double)))==NULL)
      {printf ("Error in vector allocation elon");
        exit(1);}
  if ((slon=(double*)malloc(npaths*sizeof(double)))==NULL)
      {printf ("Error in vector allocation slon");
        exit(1);}
  if ((coslat=(double*)malloc(npaths*sizeof(double)))==NULL)
      {printf ("Error in vector allocation coslat");
        exit(1);}
 
  /* allocate memory for stations*/
  if ((stalat=(float*)malloc(nsta*sizeof(float)))==NULL)
      {printf ("Error in vector allocation stalat");
        exit(1);}
  if ((stalon=(float*)malloc(nsta*sizeof(float)))==NULL)
      {printf ("Error in vector allocation stalon");
        exit(1);}
  /* end allocate memory for paths and stations*/

  for(ista=0;ista<nsta;ista++) {
    fscanf(indat,"%s %f %f ",nomsta[ista],&stalat[ista],&stalon[ista]);
    if(stalon[ista]<0.) stalon[ista]=stalon[ista]+360.;
    if((stalon[ista]<0.) || (stalon[ista]>360.)) {
      printf ("\n Pb with lon sta \n ");
      exit(1);
    }
  }

  printf(" \n Debut initialisation\n");

  for (ii=0;ii<nnodes;ii++) { 
    qualit[ii]=0;
    nopathpernode[ii]=0;
    length[ii]=0;
    qualitaz[ii][0]= '0';        
    qualitaz[ii][1]= '0';        
    qualitaz[ii][2]= '0';        
    qualitaz[ii][3]= '0';        
    qualitaz[ii][4]= '0';        
  } 
  printf(" Fin initialisation\n ");
/* 
   loops on the path, then on each point of the path, then 
   localisation of the node.
   Within the loops, we compute for each node 
   a quality criterion which depends on the length and azimuth
   of each path in each cell. We need to know :
   - the number of paths per node (inpn)
   - for each node, the length and azimuth of each paths crossing
     that nodes.
*/

  intrala=0;
  for (ipath=0;ipath<npaths;ipath++) {
    for(ik=0;ik<nnodes;ik++) {
      itest[ik]=0;
      azcell[ik]=0;
    }
    ijmem=-5;
    fgets(ligne,400,indat);
    sscanf(ligne,"%s %s ",bids,nomdatray);
    strcpy(nombid,bids);
    ilen=strlen(nombid);
    for(ii=0;ii<ilen;ii++) {
      comp=nombid[ii];
      comp2='.';
      if(comp==comp2) ichar=ii;
    }
    strcpy(nombid+ichar,"");

    /* read the event lat lon */
    fgets(ligne,400,indat);
    sscanf(ligne,"%f %f ",&ela,&elo);

    /* read and ignore the model lines */
    fgets(ligne,400,indat);
    /* read and ignore the error lines */
    fgets(ligne,400,indat);

    elon[ipath]=(double)elo;
    coelat[ipath]=(double)ela;
    if(elon[ipath]<0.) elon[ipath]=elon[ipath]+360.;
    if((elon[ipath]<0.) || (elon[ipath]>360.)) {
      printf ("\n Pb with lon path: elon= %f path %d \n ",elon[ipath],ipath);
      exit(1);
    }

    icompt=0;
    for (ista=0;ista<nsta;ista++) {
      strcpy(nombid2,nomsta[ista]);
      if(strcmp(nombid,nombid2)==0) {
        coslat[ipath]=(double)stalat[ista];
        slon[ipath]=(double)stalon[ista];
        icompt++;
      }
    }
 
    if(icompt!=1) {
      printf ("\n Pb finding station : icompt= %d \n ",icompt);
      printf("\n  ipath %d nomdatray %s nombid %s \n",ipath,nomdatray,nombid);
      exit(1);
    }
 
    coelat[ipath]=atan(0.993277*tan(coelat[ipath]*degtorad));
    coelat[ipath]=pis2-coelat[ipath];
    coslat[ipath]=atan(0.993277*tan(coslat[ipath]*degtorad));
    coslat[ipath]=pis2-coslat[ipath];
    elon[ipath]=elon[ipath]*degtorad;
    slon[ipath]=slon[ipath]*degtorad;
    delta=del(coelat[ipath],elon[ipath],coslat[ipath],slon[ipath]);
    azes=azimuth(pis2-coelat[ipath],elon[ipath],pis2-coslat[ipath],slon[ipath],delta);
    nt=(int)(8.*delta*radtodeg+1.);
    dpas=delta/(double)nt;
    nt++;
    for(jj=0;jj<nt;jj++) { 
      dst=dpas*(double)(jj);
      if (jj == 0) {
	colati=coelat[ipath];
        loni=elon[ipath];
      }
      else if (jj == nt-1) {
        colati=coslat[ipath];
        loni=slon[ipath];
      }
      else
        coord(coelat[ipath],elon[ipath],dst,azes,&colati,&loni);
           
      xd[0]=(pis2-colati)*radtodeg; 
      xd[1]=loni*radtodeg; 
      findnodesph_(&intrala,xd,&node,&iwalk,&xlon,&xlat,&dimxd);
      intrala=1;
      ij=node-1;
 
      if (itest[ij]==0) {
        nopathpernode[ij]++;
        if (ijmem>=0) {
          if (itest[ijmem]==0) {
            printf("\n itest[%d]=0\n",ijmem+1);
            exit(1);
	  } 
          azcell[ijmem]=azcell[ijmem] / (double)itest[ijmem];
          if (azcell[ijmem]>= pi || azcell[ijmem]< 0) {
            printf("\n CAS 1 azcell[%d] = %f ",ijmem+1,azcell[ijmem]);
            printf("\n itest[] =  %d\n",itest[ijmem]);
            printf("\n JJ =  %d\n",jj);
            exit(1);
	  }   
          if(azcell[ijmem]*radtodeg >= 0 && azcell[ijmem]*radtodeg <=36) 
            qualitaz[ijmem][0]= '1';
          else if(azcell[ijmem]*radtodeg > 36 && azcell[ijmem]*radtodeg <=72)
            qualitaz[ijmem][1]= '1';
          else if(azcell[ijmem]*radtodeg > 72 && azcell[ijmem]*radtodeg <=108)
            qualitaz[ijmem][2]= '1';
          else if(azcell[ijmem]*radtodeg > 108 && azcell[ijmem]*radtodeg <=144)
            qualitaz[ijmem][3]= '1';
          else if(azcell[ijmem]*radtodeg > 144 && azcell[ijmem]*radtodeg <180)
            qualitaz[ijmem][4]= '1';
          else { 
            printf("\n NODES %d azcell[%d]= %f \n",ijmem+1,ijmem,azcell[ijmem]);
            exit(1);
	  }
 
          qualit[ijmem]= qualit[ijmem] + length[ijmem];
        } 
      }
      itest[ij]++;
      length[ij]=length[ij]+dpas;

      if(jj==0)
        azim=azes;  
      else if (jj!=nt-1) {
        azim=azimuth(pis2-colati,loni,pis2-coslat[ipath],slon[ipath],delta-dst);
        if (azim >=2*pi || azim< 0) {
	  printf("\n azim %f ij+1 %d",azim,ij+1);
          exit(1);
	}
        if (azim >= pi) azim=azim-pi;
        if (azim < 0 || azim >= pi) {
	  printf("\n NODES %d  PATH %d AZIM %f ",ij+1,ipath,azim);
          printf("\n JJ %d  NT %d ",jj,nt);
	}
        azcell[ij]=azcell[ij] + azim;
        ijmem=ij;
        if(ij>64800) {
          printf("\n IJ trop grand %d ",ij);
          exit(1);
	}
        azimem=azim;
      } 
      else if (jj==nt-1) {
        azim=azimem;
        azcell[ij]=azcell[ij] + azim;
        if (itest[ij]==0) {
          printf("\n itest[%d]=0\n",ij+1);
          exit(1);
	}
        azcell[ij]=azcell[ij] / (double)itest[ij];
        if (azcell[ij]>= pi || azcell[ij]< 0) {
          printf("\n CAS 2 azcell[%d] = %f ",ij+1,azcell[ij]);
          printf("\n itest[] =  %d",itest[ij]);
          exit(1);
	}
        if(azcell[ij]*radtodeg >= 0 && azcell[ij]*radtodeg <=36)
          qualitaz[ij][0]= '1';
        else if(azcell[ij]*radtodeg > 36 && azcell[ij]*radtodeg <=72)
          qualitaz[ij][1]= '1';
        else if(azcell[ij]*radtodeg > 72 && azcell[ij]*radtodeg <=108)
          qualitaz[ij][2]= '1';
        else if(azcell[ij]*radtodeg > 108 && azcell[ij]*radtodeg <=144)
          qualitaz[ij][3]= '1';
        else if(azcell[ij]*radtodeg > 144 && azcell[ij]*radtodeg <180)
          qualitaz[ij][4]= '1';
        else {
  	  printf("\n NODES %d azcell[%d]= %f \n",ij+1,ij,azcell[ij]);
          exit(1);
	}
 
        qualit[ij]= qualit[ij] + length[ij];
      }
    }
  }        
  fclose(indat);

/* Write in the output file :lat, lon node, quality criterion  */

  printf("\n Output file \n");
  if(!(output=fopen("NODESQUAL","w"))){
    printf ("\n forcbrut: pb opening output file NODESQUAL");
    exit(1);
  }
  fprintf (output," %d %d %d",nnodes,irun,nbid);  
  for(jj=0;jj<nnodes;jj++) {
    lonnode[jj]=lonnode[jj]*radtodeg;
    colatnode[jj]=pis2-colatnode[jj];
    colatnode[jj]=radtodeg*atan(tan(colatnode[jj])/0.993277);
    fprintf (output,"\n%f %f %f ",lonnode[jj],colatnode[jj],qualit[jj]);  
    for(ii=0;ii<5;ii++) 
      fprintf (output,"%c",qualitaz[jj][ii]); 
  }
  fclose(output);
}
/*######distance between 2 points on great circle######*/
/*The argument late,lone,lats,lons are supposed to be given
  in radians. late and lats are colatitude (colatitude)*/
double del(late,lone,lats,lons) 
double late,lone,lats,lons;
{
  double resu,ca,cb,sa,sb;
 
 
  ca=cos(late);
  cb=cos(lats);
  sa=sin(late);
  sb=sin(lats);
  resu=ca*cb+(sa*sb)*cos(lone-lons);
 
  if(resu>1. || resu<-1.) 
    printf("Error forcbrut, subroutine del");
  resu=acos(resu);
  return (resu);
}
 
/* ##########great cicle azimuth########################## */
/*The argument late,lone,lats,lons,delt are supposed to be given
  in radians. late and lats are latitude (not colatitude)*/
double azimuth(late,lone,lats,lons,delt)
double late,lone,lats,lons,delt;
{
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
return (azim);
}
 
/* #################################### */
/* computation of lati and loni for a point i located on the great circle
   at a distance delt of the epicenter(late,lone). The argument late,lone,
   delt,azes are supposed to be in radians and the results lati(=colatitude)
   are in RADIANS.late and lati are colatitude (not latitude)*/

void coord(late,lone,delt,azes,lati,loni)
double late,lone,delt,azes;
double *lati,*loni;
{
double ca,cb,cc,sa,sb,sc,cosA,C,cosC;

 /*calcul de la latitude */

  cosA=cos(azes);
  cb=cos(late);
  sb=sin(late);
  cc=cos(delt);
  sc=sin(delt);

  ca=cb*cc+sb*sc*cosA;
  if(ca>1. || ca <-1.)
      {printf("Error 1 forcbrut, routine coord");
       /*exit(1);*/}
     *lati=acos(ca);

 /*calcul de la longitude*/
 /*calcul de la longitude*/

   sa=sin(*lati);

   cosC=(cc-ca*cb)/(sa*sb);
/*if(cosC<-1. || cosC>1. )
      {printf("Error 2 forcbrut, routine coord CosC = %E\n ",cosC);
         exit(1);}*/
   if (cosC>1.) cosC=1.;
   if (cosC<-1.) cosC=-1.;
   C=acos(cosC);
   if(azes>=0 && azes<=pi) *loni=C+lone;
   if(azes>pi && azes<2*pi) *loni=lone-C;
  }
/*#### memory allocation for matrice (double) #######*/
double** allotabd(l,c)
int l,c; 
{
double **mat;
int i;
if ((mat=(double**)malloc(l*sizeof(double*)))==NULL)
   {printf ("Error in matrix allocation 1");
    exit(1);}
    for (i=0;i<l;i++)
     {   
     if ((mat[i]=(double*)malloc(c*sizeof(double)))==NULL)
       {printf ("Error in matrix allocation 2");
        exit(1);}
     }   
return(mat);
}
/*#### memory allocation for matrice (int) #######*/
int** allotabi(l,c)
int l,c; 
{
int **mat;
int i;   
if ((mat=(int**)malloc(l*sizeof(int*)))==NULL)
   {printf ("Error in matrix allocation for integer 1");
    exit(1);}
    for (i=0;i<l;i++)
     {   
     if ((mat[i]=(int*)malloc(c*sizeof(int)))==NULL)
       {printf ("Error in matrix allocation for integer 2");
        exit(1);} 
     }   
return(mat);
}
