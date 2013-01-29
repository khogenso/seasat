
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
  double a, b, c;

  if (argc != 4) {printf("Usage: %s c b a\n",argv[0]); exit(1); }
  
  
  a = atof(argv[3]);
  b = atof(argv[2]);
  c = atof(argv[1]);
  
  printf("Making file for Y = %lf x^2 + %lf x + %lf\n",a,b,c); 

  for (i=0; i<LEN; i++) dop[i] = a*i*i + b*i + c;
  
  printf("Writing out caclualted dop file\n");  
  fp = fopen("dop.calc","w");
  for (i=0; i<LEN; i++) fprintf(fp,"%i %lf\n",i,dop[i]);
  fclose(fp);

}
