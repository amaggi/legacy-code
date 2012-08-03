#include"header.h"

int main(int argc, char *argv[]){

  int parsed, n_read, n_written;
  char outfile[FILENAME_MAX];

  parsed=parser(argc, argv);
  if(!parsed){
    return 1;
  }

  /* tell us what you think you are doing */
  print_parsed_info();

  /* open input and output files */
  fpin=fopen(cat_filename,"r");
  if(fpin==NULL){
    fprintf(stderr,"Cannot open %s for reading.  Exiting.\n",cat_filename);
    return 1;
  }
  sprintf(outfile,"%s.txt",basename);
  fpout=fopen(outfile,"w");
  if(fpout==NULL){
    fprintf(stderr,"Cannot open %s for writing.  Exiting.\n",outfile);
    return 1;
  }

  /* sort out headers */
  print_headers(argc, argv);

  n_read=0;
  n_written=0;
  while(!feof(fpin)){
    read_event();
    n_read++;
    if(check_criteria()){
      write_event();
      n_written++;
    }
  }

  /* close files */
  fclose(fpin);
  fclose(fpout);

  /* report */
  printf("%d / %d events matched criteria\n", n_written, n_read);

  return 0;
}
