/*
* $Id: respc.c,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
*
* $Log: respc.c,v $
* Revision 1.1.1.1  2002/07/12 11:15:20  maggi
*
*
* Revision 1.1  2002/05/23 10:28:52  maggi
* Initial revision
*
*/
/* A simple C function calling evresp.c and ready callable by a FORTRAN    */

#include <stdio.h>

#define NOERROR    0
#define NOSTATION -1
#define NOCHANNEL -2
#define NODATIME  -3
#define MAXNOUT    4

respc_(fr,sta,cha,units,file,opt,datime,prl,pim,stal,chal,unitsl,filel,optl,datimel)

char *sta,*cha,*units,*file,*opt,*datime;
int stal,chal,unitsl,filel,optl,datimel;
float *fr,*prl,*pim;
{
    int flag,i;
    float out[MAXNOUT*2];
    int nout;

/*  set up the character variables */

         for (i=0; (sta[i] != ' ') && (i < (stal-1)); ++i);
         sta[i]='\0';

         for (i=0; (cha[i] != ' ') && (i < (chal-1)); ++i);
         cha[i]='\0';

         for (i=0; (units[i] != ' ') && (i < (unitsl-1)); ++i);
         units[i]='\0';

         for (i=0; (file[i] != ' ') && (i < (filel-1)); ++i);
         file[i]='\0';

         for (i=0; (opt[i] != ' ') && (i < (optl-1)); ++i);
         opt[i]='\0';

         for (i=0; (datime[i] != ' ') && (i < (datimel-1)); ++i);
         datime[i]='\0';


    nout = 0;

        	flag = evresp_(sta,cha,datime,fr,out,&nout,units,file,opt);

		if( flag < NOERROR){
			if(opt!=NULL & flag==NOSTATION)
				printf("(station not found)\n");
			else if(opt!=NULL & flag==NOCHANNEL)
				printf("(channel not found)\n");
			else if(opt!=NULL & flag==NODATIME)
				printf("(datime not found)\n");
			else if(opt)
				printf("\n");
			exit(2);
		}


                if( nout > 1){
                        printf("(multiple entries.)\n");
                        exit(2);
                }

     *prl = out[0];
     *pim = out[1];

}
