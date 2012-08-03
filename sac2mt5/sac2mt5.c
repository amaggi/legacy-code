/*
SAC2MT5
	This program prepares sac data for input into MT5.  Input data is
	as follows:

		* Alpha-formatted sac files, with the headers corresponding
		to the actual data in the files (eg: if a broad band data 
		file has had its broad band response deconvolved and has 
		been reconvolved with a long period response, the component
		header element should read LH?, not BB?).

		* Response files for ALL the data files.  These responses
		should correspond to the data: if reconvolution has taken 
		place, then the response file should correspond to the 
		reconvolved time series.  If a response is common to more 
		than one data file, then only one copy of the response file 
		is necessary.  The program can deal with response files which 
		contain different responses for the same intstrument at 
		different times (this happens when extracting responses from 
		SEED files which contain data for more than one event). 

		* A .txt file which contains information on the event.  It
		must contain a single line in the -jdy format (any number of 
		BLANK header lines are OK).  It is particularly important that 
		the depth and origin times in the file are self consistent.
		The file should be renamed:
			input.txt

		* A file called
			sac2mt5_input
		which contains a list of filenames, of the sac data files and 
		their corresponding response files.  For example:
			sac_file1	resp_file1
			sac_file2	resp_file2
			.........       ..........
			sac_filen	resp_filen
		Note that names of response files may be repeated if they 
		correspond to more than one data file. 

	Output:
		* A file called 
			new_output
		which is formatted according to the .dsn format explained 
		in the MT5 docuentation and repeated in the routine
		write_dsn.  This file is ready to be transformed into a 
		DOS format and fed into MT5INT.

Alessia Maggi	Completed March 1999
*/

#include "sac2mt5.h" 	/* header file containing include files, 
			structure definitions and function declarations */

int main()
{

Pz poze[MAXRESP];	/* array of pole zero structures */
Sac sac_data;		/* a sac structure */
Event quake;		/* and event structure */
Npair names[MAXSTA];
int n_resp, pz_index,n_comp, i;
FILE *respfile, *sacfile, *outfile, *infile, *txtfile, *listfile;

/* 	read the .txt file 	*/
txtfile = fopen("input.txt","r");
if(txtfile==NULL){
	fprintf(stderr,"Cannot open the txt-file input.txt: aborting.\n");
	return NOTOK;
}
read_txt(txtfile,&quake);
fclose(txtfile);

/* 	open the output file 	*/
outfile = fopen("new_output","w");

/* 	write the header line of the output file 	*/
write_header(outfile,&quake);

/* 	find out how many data files we have and their names 	*/
listfile = fopen("sac2mt5_input","r");
n_comp = find_pairs(listfile,names);
fclose(listfile);


/* 	for each data file 	*/
for(i=0;i<n_comp;i++){

	/* open a pair of data file and response file */

	respfile = fopen(names[i].resp,"r");
	sacfile = fopen(names[i].data,"r");

	fprintf(stderr,"%i: Working on %s <> %s\n",i+1,names[i].data, names[i].resp);

	if(respfile != NULL && sacfile != NULL) {




	/* get the information out of the data file */
	if(read_sac(&sac_data, sacfile)==NOTOK){
		fprintf(stderr,"Error in reading the sac file: exiting.\n");
		fclose(sacfile);
		return NOTOK;
	}
	fclose(sacfile);




	/* get the information out of the response files */
	n_resp = resp2pz(poze, respfile);
	if(n_resp == NOTOK){
		fprintf(stderr,"Error in reading the response file: exiting.\n");
		fclose(respfile);
		return NOTOK;
	}
	fclose(respfile);

	/* decide which poles and zeros are applicable here */
	pz_index = which_pz(&sac_data, poze, n_resp);
	if(pz_index == NOTOK){
		fprintf(stderr,"Error in finding the correct poles and zeros :exiting.\n");
		return NOTOK;
	}

	/* write the relevant bit of the output file in dsn format*/
	if (write_dsn(outfile, &sac_data, &poze[pz_index])==NOTOK){
		fprintf(stderr,"Error in writing the output file: exiting.\n");
		fclose(outfile);
		return NOTOK;
	}

	}	
	else{	/* could not open one or both of the response and sac files */ 
		fprintf(stderr,"Couldn't open one of %s and %s\n",names[i].data, names[i].resp);
		fprintf(stderr,"Continuing with next file...\n");
	}
/* end for */
}

fclose(outfile);

}
