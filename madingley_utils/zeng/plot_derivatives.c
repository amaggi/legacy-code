#include<stdio.h>
#define LINE_MAX 120
#define MAX_LAYERS 100
#define HALF_MAX_LAYERS 50

typedef struct LAyer{
  float z, rho, vp, vs;
} Layer;

int main(int argc, char *argv[])
{
  char line[LINE_MAX];
  FILE *fpmodel, *fpin, *fpout;

  int n_layers, n_periods, i, j, dummy;
  float h, period, phase_vel, group_vel;
  float dcdvs, dcdvp, dcdrho;
  Layer model[MAX_LAYERS];

  // If number of arguments is incorrect print usage
  if(argc < 3){
    printf("Usage: plot_derivatives inputfile outputfile\n");
    printf("where inputfile is the output of dispr, and outputfile is to be read by the gmt script plot_derivatives.gmt\n");
    return;
  }

  // Tell the user which files you are reading
  printf("Reading model from dispr.in, derivatives from %s and writing to %s\n", argv[1], argv[2]);
  fpmodel=fopen("dispr.in", "r");
  fpin=fopen(argv[1], "r");
  fpout=fopen(argv[2], "w");

  // Read the model file and translate layer indices into depths
  fscanf(fpmodel, "%d", &n_layers);
  if(n_layers > HALF_MAX_LAYERS){
    printf("Increasce MAX_LAYERS and HALF_MAX_LAYERS\n");
    return;
  }
  for(i=0;i<n_layers;i++){
    // The upper depth of this layer is equal to the bottom depth 
    // of the layer above.
    if(i==0)
      model[0].z=0;
    else
      model[2*i].z=model[2*i-1].z;
    // Read the values and assign them to the upper depth of this layer
    fscanf(fpmodel, "%f %f %f %f", &h, &model[2*i].rho, &model[2*i].vp, &model[2*i].vs);
    // Compute the lower depth of this layer
    model[2*i+1].z=model[2*i].z+h;
    // Assign the same model values to this lower depth
    model[2*i+1].rho=model[2*i].rho;
    model[2*i+1].vp=model[2*i].vp;
    model[2*i+1].vs=model[2*i].vs;
  }

  // Read the number of periods
  fscanf(fpmodel,"%d", &n_periods);
  fclose(fpmodel);
  // Write the model out to file
  fpmodel=fopen("model.out", "w");
  for(i=0;i<n_layers;i++){
    fprintf(fpmodel,"%.2f %.3f %.3f %.3f\n", -1*model[2*i].z, model[2*i].rho, model[2*i].vp, model[2*i].vs);
    fprintf(fpmodel,"%.2f %.3f %.3f %.3f\n", -1*model[2*i+1].z, model[2*i+1].rho, model[2*i+1].vp, model[2*i+1].vs);
  }
  fclose(fpmodel);

  // Read the derivative information from the input file and write 
  // to output file
  for(i=0;i<n_periods;i++){
    // Read two blank lines
    fgets(line, LINE_MAX, fpin);
    //printf("1) %s", line);
    fgets(line, LINE_MAX, fpin);
    //printf("2) %s", line);
    // Read period, phase velocity, group velocity
    fscanf(fpin,"%f %f %f", &period, &phase_vel, &group_vel);
    // Print a separation character into the output file
    //fprintf(fpout, "> -Z%.3f T=%.3f C=%.3f U=%.3f\n", period, period, phase_vel, group_vel);
    // Read another blank lines
    fgets(line, LINE_MAX, fpin);
    //printf("3) %s", line);
    fgets(line, LINE_MAX, fpin);
    //printf("4) %s", line);
    fgets(line, LINE_MAX, fpin);
    //printf("5) %s", line);
    //Read the derivatives and write to file
    for(j=0;j<n_layers;j++){
      fscanf(fpin, "%d %f %f %f",&dummy, &dcdvs, &dcdvp, &dcdrho);
      //printf("%d %f %f %f\n", dummy, dcdvs, dcdvp, dcdrho);
      fprintf(fpout,"%.2f %.5f %.5f %.5f %.2f \n", -1*model[2*j].z, dcdrho, dcdvp, dcdvs, period);
      fprintf(fpout,"%.2f %.5f %.5f %.5f %.2f \n", -1*model[2*j+1].z, dcdrho, dcdvp, dcdvs, period);
    }
    fgets(line, LINE_MAX, fpin);
    //printf("6) %s", line);
  }

  fclose(fpin);

  //fprintf(fpout, ">\n");
  fclose(fpout);

}

