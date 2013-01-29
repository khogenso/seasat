
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LEN 6840

main (int argc, char *argv[])
{
  FILE *fp;
  int line[LEN];
  double dop[LEN];
  int i;
  
  printf("Reading in dop file\n");
  fp = fopen("dop.out","r");
  for (i=0; i<LEN; i++) fscanf(fp,"%i %lf\n",&line[i],&dop[i]);
  fclose(fp);
  
  for (i=1;i<LEN;i++) {
    if (fabs(dop[i]-dop[i-1]) > 0.7) {
      if ((dop[i]-dop[i-1])>0.0) dop[i] = dop[i]-1.0;
      else dop[i] = dop[i]+1.0;
    }
  }
 
  printf("Writing out new dop file\n");  
  fp = fopen("dop.new","w");
  for (i=0; i<LEN; i++) fprintf(fp,"%i %lf\n",line[i],dop[i]);
  fclose(fp);
}

  
