/*
SAC2MT5
	which_pz
		This routine chooses which out of a set of pole_zero structures
		contains the correct response for the time of the event.  It
		assumes that the responses have been read out of the file in 
		CRONOLOGICAL order.  The routine is passed a pointer to a
		Sac structure, an array of Pz structures and the number of
		Pz structures which are going to be looked at.  It returns 
		the index [0 -> (n-1)] of the array of Pz structures which 
		contains the correct pole zero response, or NOTOK if it 
		encounters a problem.

Alessia Maggi	Completed March 1999
*/


#include "sac2mt5.h"

int which_pz(Sac *data, Pz poze[], int n_pz)
{

int i;

switch (n_pz){
	case 0 : 	/* No responses are present */
		fprintf(stderr,"There are no responses!!\n");
		return NOTOK;
		break;
	case 1 :  	/* There is only one response: either it is OK or
			it is not */
		if(data->nzyear < poze[0].year && data->nzjday < poze[0].jday){
			fprintf(stderr,"Response is too late!\n");
			return NOTOK;
		}
		else
			return 0;
		break;
	default :	/* check if responses are too late */ 
		if(data->nzyear < poze[0].year && data->nzjday < poze[0].jday){
			fprintf(stderr,"Response is too late!\n");
			return NOTOK;
		}
			/* check last response first */ 
		else if(data->nzyear >= poze[n_pz-1].year && data->nzjday >= poze[n_pz-1].jday)
		return n_pz-1;
			/* the best pole zero set is the last one */
		else 	/* check all other responses */
			for(i=0;i<n_pz-1;i++){
				if(data->nzyear >= poze[i].year && data->nzjday >= poze[i].jday && data->nzyear < poze[i+1].year && data->nzjday < poze[i+1].jday)
					return i;
			}
		break;
	
}

/* If you get here something really funny has happened */
fprintf(stderr,"Trouble in which_pz.c\n");
return NOTOK;
}
