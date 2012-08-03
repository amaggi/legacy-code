#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>

/* delete a specified amount of node according  
   the quality criterion for azimuth read in the file
   NODESQUAL. Quality of the node is sorted according to
   the following criteria :
         case 11111 : class[ii]=1;

         case 11110 :
         case 11101 :
         case 11011 :
         case 10111 :
         case 1111 : class[ii]=2;     01111

         case 10101 :
         case 1011 :                  01011
         case 11010 :
         case 1101 :                  01101
         case 10110 : class[ii]=3;

         case 111 :                   00111
         case 10011 :
         case 11001 :
         case 11100 :
         case 1110 : class[ii]=4;     01110

         case 10100 :
         case 1010 :                  01010
         case 101 :                   00101
         case 10010 :
         case 1001 : class[ii]=5;     01001

         case 11000 :
         case 1100 :                   01100
         case 110 :                    00110
         case 11 :                     00011
         case 10001 : class[ii]=6;

         case 10000 :
         case 1000 :                   01000
         case 100 :                    00100
         case 10 :                     00010
         case 1 : class[ii]=7;         00001

         case 0 : class[ii]=8;         00000

*/

/* #### Functions and void #### */

void randelnode();
float ran3();
/* ####  MAIN PGM   ################### */
main()
{
char ligne[100],nomout[100],no[10];
int *class,*qualitaz,*statnode;
float *lonnode,*latnode,*qualit;
float prop;
int nnodes,ii,total,totdel,delnode,statusnode,bid,iter;
int iclass,idelnode,nnodesnew,nnodesinit;
int compt[8],idtot[8];
FILE *inodes;
FILE *onodes;

  for (ii=0;ii<8;ii++)
    {   idtot[ii]=0;
        compt[ii]=0; }
 
  printf ("\n Proportion of nodes to be deleted : \n");
  gets(ligne);
  sscanf(ligne," %f",&prop);
  
  if(!(inodes=fopen("NODESQUAL.work","r")))
    {printf ("\n pb opening file NODESQUAL.work");
     exit(1);}
  fgets(ligne,100,inodes);
  sscanf(ligne,"%d %d %d",&nnodes,&bid,&nnodesinit);
/*  delnode=(int)nnodesinit*prop/100.; */
  delnode=(int)nnodes*prop/100.;
  
  if ((latnode=(float*)malloc(nnodes*sizeof(float)))==NULL)
    {printf ("Error in vector allocation latnode");
     exit(1);}
  if ((lonnode=(float*)malloc(nnodes*sizeof(float)))==NULL)
    {printf ("Error in vector allocation lonnode");
     exit(1);}
  if ((qualit=(float*)malloc(nnodes*sizeof(float)))==NULL)
    {printf ("Error in vector allocation qualit");
     exit(1);}
  if ((qualitaz=(int*)malloc(nnodes*sizeof(int)))==NULL)
    {printf ("Error in vector allocation qualitaz");
     exit(1);}
  if ((class=(int*)malloc(nnodes*sizeof(int)))==NULL)
    {printf ("Error in vector allocation class");
     exit(1);}
  if ((statnode=(int*)malloc(nnodes*sizeof(int)))==NULL)
    {printf ("Error in vector allocation statnode");
     exit(1);}

  for (ii=0;ii<nnodes;ii++)
    {
    statnode[ii]=0;
    fgets(ligne,100,inodes);
    sscanf(ligne,"%f %f %f %d",&lonnode[ii],&latnode[ii],&qualit[ii],&qualitaz[ii]);
/*  printf("II= %d %f %f %f %d\n",ii,lonnode[ii],latnode[ii],qualit[ii],qualitaz[ii]); */
    switch(qualitaz[ii])
        {
         case 11111 : class[ii]=0;
                      compt[0]++;
                      break;
         case 11110 :
         case 11101 :
         case 11011 :
         case 10111 :
         case 1111 : class[ii]=1;
                      compt[1]++;
                      break;
         case 10101 :
         case 1011 :
         case 11010 :
         case 1101 :
         case 10110 : class[ii]=2;
                      compt[2]++;
                      break;
         case 111 :
         case 10011 :
         case 11001 :
         case 11100 :
         case 1110 : class[ii]=3;
                      compt[3]++;
                      break;
         case 10100 :
         case 1010 :
         case 101 :
         case 10010 :
         case 1001 : class[ii]=4;
                      compt[4]++;
                      break;
         case 11000 :
         case 1100 :
         case 110 :
         case 11 :
         case 10001 : class[ii]=5;
                      compt[5]++;
                      break;
         case 10000 :
         case 1000 :
         case 100 :
         case 10 :
         case 1 : class[ii]=6;
                      compt[6]++;
                      break;
     
         case 0 : class[ii]=7;
                      compt[7]++;
                      break;
         default  : printf("QUALITAZ[%d]= %d\n",ii,qualitaz[ii]);
         }
    }
  printf (" Initial number of Nodes %d:\n",nnodes);
  printf (" Number of Nodes in :\n");
  printf (" -class 1 : %d\n",compt[0]);
  printf (" -class 2 : %d\n",compt[1]);
  printf (" -class 3 : %d\n",compt[2]);
  printf (" -class 4 : %d\n",compt[3]);
  printf (" -class 5 : %d\n",compt[4]);
  printf (" -class 6 : %d\n",compt[5]);
  printf (" -class 7 : %d\n",compt[6]);
  printf (" -class 8 : %d\n",compt[7]);
  fclose(inodes);
 total=0;
 totdel=0;
 iclass=7;
 iter=8;
 while(iclass>0)
   {
    total=total+compt[iclass];
 /* printf (" total iclass no %d  : %d \n",iclass,total); */
     if ( (delnode-compt[iclass])>=0  && (iclass==1) && (delnode!=0) )
          {
      /*   if((int)compt[iclass]/2 >= (int)nnodesinit/100) */
           if((int)compt[iclass]/2 >= (int)nnodesinit/100) 
             {
              idtot[iclass]=1;
              delnode=(int)compt[iclass]/2;
              totdel=totdel+delnode;
              randelnode(iclass,nnodes,delnode,class,statnode);
              iter=iclass+1;
             }
           else
             {
              idtot[iclass]=2;
              delnode=compt[iclass];
              totdel=totdel+delnode;
              printf("\n The %d nodes in class %d will be deleted ",delnode,iclass+1);
              printf("\n before to exit\n");
              iter=1;
             }
          }
     else if ( (delnode-compt[iclass])>=0  && (iclass!=1) && (delnode!=0) )
          {
           idtot[iclass]=2;
           totdel=totdel+compt[iclass];
           delnode=delnode-compt[iclass];
           printf("\n The %d nodes in class %d will be deleted \n",compt[iclass],iclass+1);
           if(delnode==0) iter=iclass+1;
          }

     else if ( (delnode-compt[iclass])<0  && (delnode!=0) )
          {
           idtot[iclass]=1;
           randelnode(iclass,nnodes,delnode,class,statnode);
           totdel=totdel+delnode;
           delnode=0;
           iter=iclass+1;
          }
     else if (delnode==0)
          {
           idtot[iclass]=0;
          }
   iclass--;
   }

  total=total+compt[0];
  printf ("    Total  : %d \n",total);
  if(total!=nnodes)
    {printf("Error in sorting nodes\n");
     exit(1);}

strcpy(nomout,"NODESQUAL.");
sprintf(no,"%d",iter);
strcat(nomout,no);
if(!(onodes=fopen(nomout,"w")))
    {printf ("\n pb opening file out %s \n",nomout);
     exit(1);}

nnodesnew=nnodes-totdel;
fprintf(onodes,"%d %d %d \n",nnodesnew,bid,nnodesinit);
statusnode=0;
for(ii=0;ii<nnodes;ii++)
  {
  iclass=class[ii];
  if (idtot[iclass]==2)
         statusnode=1;
  if (idtot[iclass]==1)
         statusnode=statnode[ii];
  if (idtot[iclass]==0)
         statusnode=0;
  if(statusnode==0)
    fprintf(onodes,"%f %f %f %d\n",lonnode[ii],latnode[ii],qualit[ii],qualitaz[ii]);
  }
fclose(onodes);
}
/*############randelnode ########## */
void randelnode(iclass,nnodes,idelclass,class,statnode)
int iclass,nnodes,idelclass;
int *class, *statnode;
{
int *tab;
int ii,jj,ij,kk,all;
long idum;

if ((tab=(int*)malloc(nnodes*sizeof(int)))==NULL)
    {printf ("Error in vector allocation tab");
     exit(1);}
jj=0;
for (ii=0;ii<nnodes;ii++)
  {
   if (class[ii]==iclass)
       {tab[jj]=ii;
        jj++;}
  }
kk=0;
   printf("\n %d nodes will be deleted  among %d in class %d\n",idelclass,jj,iclass+1);
idum=78061433;
/* keep in mind that (int)23.99=23
   and (int)23.1=23*/
for (ij=0;ij<idelclass;ij++)
  {
/* all=(int)((jj+1)*ran3(&idum));  jj+1 faux, on met jj a cause du jj++ plus haut*/
   all=(int)(jj*ran3(&idum));
   if(ij==10000) idum=(long)(-all*100000000000); /* after 10000 call the seed is reinitialised to
                                                    start with a new sequence */
   kk=tab[all];
/* printf("\n ALL = %d KK = %d IJ = %d \n",all,kk,ij); */
   if(statnode[kk]==1)ij--;
   statnode[kk]=1;
  }
}
/*############ran3 generate a random number between 0 and 1########## */
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
float ran3(long *idum)
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/* (C) Copr. 1986-92 Numerical Recipes Software jn>#(3#(11,1&G2v3#)K. */

