#include<stdio.h>
#include<stdlib.h>
#include<string.h>

char scratch_string[512];

void copy_to_scratch_and_terminate(const char *string, int len){

  strncpy(scratch_string,string,len);
  scratch_string[len]='\0';

}

int fixed_len_atoi(const char *string, int len){

  copy_to_scratch_and_terminate(string,len);
  return atoi(scratch_string);
 
}

double fixed_len_atof(const char *string, int len){

  copy_to_scratch_and_terminate(string,len);
  return atof(scratch_string);

}

