/******************************************************************************
NAME: fix headers - cleans up a messy header file for SEASAT processing

SYNOPSIS: fix_headers <infile> <outfile>

DESCRIPTION:
	<infile> is the input header file
	<outfile> if the output header file after cleansing


EXTERNAL ASSOCIATES:
    NAME:               USAGE:
    ---------------------------------------------------------------

FILE REFERENCES:
    NAME:               USAGE:
    ---------------------------------------------------------------

PROGRAM HISTORY:
    VERS:   DATE:  AUTHOR:      PURPOSE:
    ---------------------------------------------------------------
    1.0	    10/12   T. Logan     Seasat Proof of Concept Project - ASF
    
HARDWARE/SOFTWARE LIMITATIONS:

ALGORITHM DESCRIPTION:

The order of columns in the header file is:

    major_cnt
    major_sync_loc,
    station_code,		set to a fixed value
    lsd_year,			set to 8
    day_of_year,		set to moving median value
    msec,			set to smoothly increasing values
    clock_drift,		?
    no_scan_indicator_bit,	
    bits_per_sample,		set to 5
    mfr_lock_bit,		
    prf_rate_code,		set to 4
    delay,			set to moving median value
    scu_bit,
    sdf_bit,
    adc_bit,
    time_gate_bit,
    local_prf_bit,
    auto_prf_bit,
    prf_lock_bit,
    local_delay_bit

ALGORITHM REFERENCES:

BUGS:

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BITS_PER_SAMPLE 5
#define PRF_RATE_CODE   4
#define LSD_YEAR        8
#define WINDOW_SIZE	100
#define RECALC_SIZE     (WINDOW_SIZE/5)

#define MAX_DOY	  400		/* day of year, but close enough!       */
#define MAX_DELAY 64    	/* one extra because we are not 1 based */
#define MAX_CLOCK_DRIFT  4097	/* 12 bit field */
#define TOLERANCE  513		/* how far off a time value can be from estimate */
#define PRI  0.60716454159

#define DISPLAY_FITS 0  /* set to 1 if you want to see each linear fit */
#define SAVE_FITS    0  /* set to 1 if you want to create the line_fits.txt file */

typedef struct {
        int       major_cnt;
	long int  major_sync_loc;
	int       lsd_year;
	int       station_code;
	long int  msec;
	int       day_of_year;
	int       clock_drift;  
	int       no_scan_indicator_bit;
	int       bits_per_sample;
	int       mfr_lock_bit;
	int       prf_rate_code;
	int       delay;
	int       scu_bit;
	int       sdf_bit;
	int       adc_bit;
	int       time_gate_bit;
	int       local_prf_bit;
	int       auto_prf_bit;
	int       prf_lock_bit;
	int       local_delay_bit;
}  SEASAT_header_ext;

int get_values(FILE *fp,SEASAT_header_ext *s);
int bitfix(double sdiff,long int *msec);
int get_median(int *hist, int size);
long int get_true_median(long int *a);
void yaxb(double x_vec[], double y_vec[],int n, double * a,double * b);

main(int argc, char *argv[])
{
  FILE *fpin,*fpout,*fpfit;
  SEASAT_header_ext **hdr;
  int i, j, k, vals;
  int icnt = 0, ocnt = 0, curr = 0, optr = WINDOW_SIZE/2;
  int fcnt = 0, bit_cnt = 0, fill_cnt = 0, line_cnt = 0;
 
  long int msec;
  double times[WINDOW_SIZE];
  double lines[WINDOW_SIZE];
  double a, b, c;
  double diff, sdiff, tmp;
  
  int offset=0, bad_cnt=0;

  if (argc!=3) {
    printf("Usage: %s <in_header_file> <out_cleaned_header_file>\n\n",argv[0]);
    printf("<in>\tName of input header file to clean");
    printf("<out>\tOutput name of cleaned header file");
    printf("\n\n");
    exit(1);
  }
  
  printf("\n\nDECODED SEASAT TIME CLEANSING PROGRAM\n\n");
  printf("\topening files...\n");
  fpin=fopen(argv[1],"r");
  if (fpin==NULL) {printf("ERROR: Unable to open input file %s\n",argv[1]); exit(1);}
  fpout=fopen(argv[2],"w");
  if (fpout==NULL) {printf("ERROR: Unable to open output file %s\n",argv[2]); exit(1);}
  if (SAVE_FITS==1) {
    fpfit=fopen("line_fits.txt","w");
    if (fpfit==NULL) {printf("ERROR: Unable to open output file line_fits.txt\n"); exit(1);}
  }
  
  
  /* Allocate the array of headers */
  hdr = (SEASAT_header_ext **) malloc(sizeof(SEASAT_header_ext *)*WINDOW_SIZE);
  if (hdr == NULL) {printf("ERROR: unable to allocate hdr\n"); exit(1);}
  for (i=0; i<WINDOW_SIZE; i++)
    hdr[i] = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));
  for (i=0; i<WINDOW_SIZE; i++)
    if (hdr[i]==NULL) {printf("Error in malloc\n"); exit(1);}

  /* Read in the first WINDOW_SIZE values and add them into the histograms
   ======================================================================*/
  for(i=0; i<WINDOW_SIZE; i++) { 
    vals = get_values(fpin,hdr[i]);
    if (vals != 20) {printf("ERROR: can't read from input file\n"); exit(1);}
    icnt++;
    times[i] = hdr[i]->msec;
    lines[i] = hdr[i]->major_cnt;
  }
  yaxb(lines,times,WINDOW_SIZE,&a,&b);

  /* Now, dump out the first WINDOW_SIZE/2 values to output file (corrected)
   ========================================================================*/
  printf("\tdumping initial lines to output file\n");
  for(i=0; i<WINDOW_SIZE/2; i++) {
  
    tmp = a*hdr[i]->major_cnt+b;
    if (fabs(tmp - hdr[i]->msec) > TOLERANCE) {
      msec = (long int) (tmp+0.5);
      printf("At %i bad value %li fixed value %li\n",hdr[i]->major_cnt,hdr[i]->msec,msec);
      hdr[i]->msec = msec;
      fcnt++;
    } else { msec = hdr[i]->msec; }
	    
    fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
      hdr[i]->major_cnt,
      hdr[i]->major_sync_loc,
      hdr[i]->station_code,
      hdr[i]->lsd_year,
      hdr[i]->day_of_year,
      msec,
      hdr[i]->clock_drift,
      hdr[i]->no_scan_indicator_bit,
      hdr[i]->bits_per_sample,
      hdr[i]->mfr_lock_bit,
      hdr[i]->prf_rate_code,
      hdr[i]->delay,
      hdr[i]->scu_bit,
      hdr[i]->sdf_bit,
      hdr[i]->adc_bit,
      hdr[i]->time_gate_bit,
      hdr[i]->local_prf_bit,
      hdr[i]->auto_prf_bit,
      hdr[i]->prf_lock_bit,
      hdr[i]->local_delay_bit);
    
    ocnt++;
  }
  
  while (vals == 20) {
    if ((icnt%10000)==0) {printf("\tcleaning line %i\n",icnt);}
  
    /* Redo the time linear regression every so often 
     -----------------------------------------------*/
    if ( (icnt%RECALC_SIZE) == 0 ) {
      // printf("Recalculating linear fit for MSEC now\n");
      for (i=0; i<WINDOW_SIZE; i++) {
        times[i] = hdr[i]->msec;
  	lines[i] = hdr[i]->major_cnt;
      }
      double old_a = a;
      double old_b = b;
      if (DISPLAY_FITS==1) printf("ICNT %i: ",icnt);
      yaxb(lines,times,WINDOW_SIZE,&a,&b);
      if (a > 0.6073 || a < 0.607 ) {  /* BAD do not use*/ 
        if (old_a > 0.6073 || old_a < 0.607) { /* last were bad too! */
          if (fabs(old_a-PRI)<fabs(a-PRI)) { /* this is worse, don't use it */
            a=old_a;
	    b=old_b; 
	    if (DISPLAY_FITS==1) printf("\tDiscarded\n");
	  } else { 
	    if (DISPLAY_FITS==1) printf("\tUsed\n"); 
	    if (SAVE_FITS==1) fprintf(fpfit,"OK: %li %lf %lf\n",icnt,a,b);  
	    offset = 0; 
	  }
	} else {
          a=old_a;
	  b=old_b; 
	  if (DISPLAY_FITS==1) printf("\tDiscarded\n");
	}
      } else { 
	if (DISPLAY_FITS==1) printf("\tUsed\n"); 
	if (SAVE_FITS==1) fprintf(fpfit,"GOOD: %li %lf %lf\n",icnt,a,b);  
	offset = 0; 
      }
    }
  
    /* read and add in the next set of values 
     ---------------------------------------*/
    vals = get_values(fpin,hdr[curr]);
    if (vals == 20) {
      icnt++;
      tmp = a*(hdr[optr]->major_cnt+offset)+b;
      diff = fabs(tmp - hdr[optr]->msec);
      sdiff = tmp - hdr[optr]->msec;
      msec = hdr[optr]->msec;
      
      if (diff > TOLERANCE) {
        if (bitfix(sdiff,&msec) == 1) {
	  printf("At %i bad value %li fixed value %li diff %lf (bit fix fill)\n",
	    hdr[optr]->major_cnt,hdr[optr]->msec,msec,sdiff);
	    bit_cnt++; fcnt++;
	} else {
          int minus1 = (optr-1+WINDOW_SIZE)%WINDOW_SIZE;
	  int plus1  = (optr+1)%WINDOW_SIZE;
          if (hdr[minus1]->msec!=hdr[optr]->msec &&    /* this is not the same as last */
	      hdr[minus1]->msec==hdr[plus1]->msec &&   /* next is the same as last     */
	      ((double)hdr[minus1]->msec-tmp)<20) {    /* within +20 of the linear trend */
 	    msec = hdr[minus1]->msec;  // fills in gaps in flat lines
	    printf("At %6i bad value %9li fixed value %9li diff %9li (gap fill %9.6lf)\n",
	      hdr[optr]->major_cnt,hdr[optr]->msec,msec,hdr[optr]->msec-msec,((double)hdr[minus1]->msec-tmp) );
	    fill_cnt++; fcnt++;
          } else {
            msec = (long int) (tmp+0.5); 
	    printf("At %i bad value %li fixed value %li diff %li (linear fill)\n",
	      hdr[optr]->major_cnt,hdr[optr]->msec,msec,hdr[optr]->msec-msec);
	    line_cnt++; fcnt++;
	  }
	}
	// hdr[optr]->msec = msec; /* fix the past points for the next fit */
      } 
      
      curr = (curr+1)%WINDOW_SIZE;
       
      fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
        hdr[optr]->major_cnt,hdr[optr]->major_sync_loc,hdr[optr]->station_code,
        hdr[optr]->lsd_year,hdr[optr]->day_of_year,msec,
        hdr[optr]->clock_drift,hdr[optr]->no_scan_indicator_bit,hdr[optr]->bits_per_sample,
        hdr[optr]->mfr_lock_bit,hdr[optr]->prf_rate_code,hdr[optr]->delay,
        hdr[optr]->scu_bit,hdr[optr]->sdf_bit,hdr[optr]->adc_bit,hdr[optr]->time_gate_bit,
        hdr[optr]->local_prf_bit,hdr[optr]->auto_prf_bit,hdr[optr]->prf_lock_bit,
        hdr[optr]->local_delay_bit);

      optr = (optr+1)%WINDOW_SIZE;
      ocnt++;
    }
  }

  /* Finally, dump out the last WINDOW_SIZE/2 values (corrected)
   =============================================================*/
  printf("\tdumping final lines to output file\n");
  
  for(i=0; i<WINDOW_SIZE/2; i++) {

    tmp = a*hdr[optr]->major_cnt+b;
    if (fabs(tmp - hdr[optr]->msec) > TOLERANCE) {
      msec = (long int) (tmp+0.5);
      printf("At %i bad value %li fixed value %li\n",hdr[optr]->major_cnt,hdr[optr]->msec,msec);
      hdr[optr]->msec = msec;
      fcnt++;
    } else { msec = hdr[optr]->msec; }

    fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
      hdr[optr]->major_cnt,
      hdr[optr]->major_sync_loc,
      hdr[optr]->station_code,
      hdr[optr]->lsd_year,
      hdr[optr]->day_of_year,
      msec,
      hdr[optr]->clock_drift,
      hdr[optr]->no_scan_indicator_bit,
      hdr[optr]->bits_per_sample,
      hdr[optr]->mfr_lock_bit,
      hdr[optr]->prf_rate_code,
      hdr[optr]->delay,
      hdr[optr]->scu_bit,
      hdr[optr]->sdf_bit,
      hdr[optr]->adc_bit,
      hdr[optr]->time_gate_bit,
      hdr[optr]->local_prf_bit,
      hdr[optr]->auto_prf_bit,
      hdr[optr]->prf_lock_bit,
      hdr[optr]->local_delay_bit);

    optr = (optr+1)%WINDOW_SIZE;
    ocnt++;
  }
  
  if (icnt != ocnt) printf("ERROR: input/output don't match; read %i wrote %i\n",icnt,ocnt);
  else {
    printf("\n\nDone with calculations - read and wrote %i lines; fixed %i time values (%f%%)\n\n",
  	icnt,fcnt, 100.0*(float)fcnt/(float)icnt);
    printf("\tTYPE      \tNumber\tChanges\tTotal\n");
    printf("\tbit fixes \t%i\t%5.2f%%\t%5.2f%%\n",bit_cnt,100*(float)bit_cnt/(float)fcnt,100*(float)bit_cnt/(float)icnt);
    printf("\tfill fixes\t%i\t%5.2f%%\t%5.2f%%\n",fill_cnt,100*(float)fill_cnt/(float)fcnt,100*(float)fill_cnt/(float)icnt);
    printf("\tline fixes\t%i\t%5.2f%%\t%5.2f%%\n\n",line_cnt,100*(float)line_cnt/(float)fcnt,100*(float)line_cnt/(float)icnt);
  }
  
  fclose(fpin);
  fclose(fpout);
  if (SAVE_FITS==1) fclose(fpfit);

  exit(0);
}

int get_values(FILE *fp,SEASAT_header_ext *s)
{
  int val;
  if (s==NULL) {printf("empty pointer passed to get_values\n"); exit(1);}
  val = fscanf(fp,"%i %li %i %i %i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    &(s->major_cnt),
    &(s->major_sync_loc),
    &(s->station_code),
    &(s->lsd_year),
    &(s->day_of_year),
    &(s->msec),
    &(s->clock_drift),
    &(s->no_scan_indicator_bit),
    &(s->bits_per_sample),
    &(s->mfr_lock_bit),
    &(s->prf_rate_code),
    &(s->delay),
    &(s->scu_bit),
    &(s->sdf_bit),
    &(s->adc_bit),
    &(s->time_gate_bit),
    &(s->local_prf_bit),
    &(s->auto_prf_bit),
    &(s->prf_lock_bit),
    &(s->local_delay_bit));

  return(val);
}

int bitfix(double sdiff,long int *msec)
{
  long int target=16;
  int n;
  for (n=4; n<34; n++) {
    target = target * 2;
    // printf("\tchecking %lf\n",sdiff);
    if (fabs(fabs(sdiff)-(double)target)<2.0) {
      // printf("Found bit fix:  original value: %li  diff: %lf  target: %li ",*msec,sdiff,target);
      if (sdiff<0) *msec=*msec-target;
      else *msec=*msec+target;
      // printf("new %li\n",*msec);
      return(1);
    }
  }
  return(0);
}

int get_median(int *hist, int size) {
  int retval = -1, max = 0, i;
  for (i=0; i<size; i++) if (hist[i]>max) {max=hist[i]; retval=i;}
  if (retval==-1) { printf("Error getting histogram median value\n"); exit(1); }
  return(retval);    
}

long int get_true_median(long int *a)
{
  long int sorted[WINDOW_SIZE];
  int i,j;
  long int t;  
  
  for(i=0;i<WINDOW_SIZE;i++) sorted[i] = a[i];
  for(i=WINDOW_SIZE-2;i>=0;i--) {  
   for(j=0;j<=i;j++) {  
    if(sorted[j]>sorted[j+1]) {  
      t=sorted[j];  
      sorted[j]=sorted[j+1];  
      sorted[j+1]=t;  
    }  
   }  
  }  

//  printf("SORTED: ");
//  for (i=0;i<WINDOW_SIZE;i++) { printf("%li ",sorted[i]); }
//  printf("\n");

  return(sorted[WINDOW_SIZE/2-1]);
}


/******************************************************************************
NAME: yaxb.c

SYNOPSIS:	yaxb(double x_vec[], double y_vec[], int n, double *a, double *b)

DESCRIPTION:	Computes a and b for y = ax + b using linear regression
		given double vectors y and x of length n.

PARAMETERS: 	x_vec   double[]     Input vector of X values
        	y_vec   double[]     Input vector of Y values
        	n   	int         Length of input vectors
        	a   	double*      Return x coefficient
        	b   	double*      Return offset factor

HISTORY:       Borowed from ASF tools and converted to double - T. Logan 10/12
	`	
ALGORITHM REF:  Cheney, Ward & D. Kincaid, Numerical Mathematics and Computing,
            2nd Editn. pp 360-362. Brooks/Cole Pub. Co., Pacific Grove, Ca.
******************************************************************************/

void yaxb(double x_vec[], double y_vec[],int n, double * a,double * b)
{
 double sum_x, sum_xx, sum_xy, sum_y;
 double d, at, bt;
 int   i, cnt;
 double res = 100.0;
 double max_val, tmp, diff;
 int    max_loc;
 
 while (res > 0.01) {
   sum_x = 0.0; sum_xx = 0.0;
   sum_y = 0.0; sum_xy = 0.0;
   cnt = 0;
   
   for (i=0; i<n; i++) {
      if (x_vec[i] != 0.0) {
        sum_x += x_vec[i]; sum_y += y_vec[i];
        sum_xx += x_vec[i] * x_vec[i];
        sum_xy += x_vec[i] * y_vec[i];
        cnt++;
      }
   }

   d =  cnt * sum_xx - sum_x*sum_x;
   at = cnt * sum_xy - sum_x*sum_y;
   bt = sum_xx * sum_y - sum_x * sum_xy;
   *a = at/d;
   *b = bt/d;
   max_val = 0.0;
   max_loc = -1;
 
   /* check the regression, throwing out the one input that has the highest error */
   for (i=0; i<n; i++) {
     if (x_vec[i] != 0.0) {
       tmp = *a * x_vec[i] + *b; 
       diff = abs(y_vec[i]-tmp);
       if (diff > max_val) { max_val = diff; max_loc = i; }
     }
   }
   
   if (max_loc != -1) {
//     printf("Culling point %i (%lf) res=%lf\n",max_loc,y_vec[max_loc],max_val);
     x_vec[max_loc]=0.0;
     y_vec[max_loc]=0.0;
     res = max_val;
   }
   else  res = 0.0;
 }
 
 if (isnan(*a) || isnan(*b)) (*a) = (*b) = 0.0;
 if (DISPLAY_FITS==1) printf("Found coefficients y = %lf x + %lf ",*a,*b);
}


