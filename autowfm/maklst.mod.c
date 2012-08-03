#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
char nombid[100];
char ligne[100];
char nomin[100];
char nomsort[100];
int npaths,i;
FILE *sortie;
FILE *entree;

main()
{
strcpy(nomin,"bid");
strcpy(nomsort,"lst.mod");
if(!(entree=fopen(nomin,"r")))
       {printf ("%s\n : pb opening file",nomin);
        exit(1);}
if(!(sortie=fopen(nomsort,"w")))
       {printf ("%s\n : pb opening file",nomsort);
        exit(1);}

fprintf (sortie, "#");
fgets(ligne,100,entree);
sscanf (ligne, "%d ", &npaths);
fprintf (sortie, "\n set NBSEI = %d", npaths);
for(i=1;i<=npaths;i++)
     {
      fgets(ligne,100,entree);
      sscanf (ligne, "%s",nombid);
      fprintf (sortie, "\n set Fic%d = %s",i,nombid);
     }
fclose (sortie);
fclose (entree);
}

