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
#define WINDOW_SIZE	400

#define MAX_DOY	  400		/* day of year, but close enough!       */
#define MAX_DELAY 64    	/* one extra because we are not 1 based */
#define MAX_CLOCK_DRIFT  4097	/* 12 bit field */
#define TOLERANCE  2		/* how far off a time value can be from estimate */
#define SHIFT_GAP  5


typedef struct {
        int      major_cnt;
	long int  major_sync_loc;
	int     lsd_year;
	int     station_code;
	long int msec;
	int      day_of_year;
	int      clock_drift;  
	int     no_scan_indicator_bit;
	int     bits_per_sample;
	int     mfr_lock_bit;
	int     prf_rate_code;
	int     delay;
	int     scu_bit;
	int     sdf_bit;
	int     adc_bit;
	int     time_gate_bit;
	int     local_prf_bit;
	int     auto_prf_bit;
	int     prf_lock_bit;
	int     local_delay_bit;
}  SEASAT_header_ext;

int get_values(FILE *fp,SEASAT_header_ext *s);
int get_median(int *hist, int size);
void yaxb(double x_vec[], double y_vec[],int n, double * a,double * b);
void yax2bxc(double x_vec[],double y_vec[],int n,double *a,double *b,double *c);

main(int argc, char *argv[])
{
  FILE *fpin,*fpout;
  FILE *fpdis;
  SEASAT_header_ext **hdr;
  int i, j, k, vals;
  int icnt = 0, ocnt = 0, curr = 0, optr = WINDOW_SIZE/2;
 
  int station_code_hist[16];    
  int doy_hist[400];
  int clock_drift_hist[MAX_CLOCK_DRIFT];
  int delay_hist[MAX_DELAY];
  int lsd_year_hist[16];
  int bits_per_sample_hist[16]; 
  int prf_rate_code_hist[16];   

  int station_code_median;
  int doy_median;
  int clock_drift_median;
  int lsd_year_median;
  int bits_per_sample_median;
  int prf_rate_code_median;
  int delay_median;
  
  long int msec;
  double times[WINDOW_SIZE];
  double lines[WINDOW_SIZE];
  double a, b, c;
  double diff, save_diff;
  
  int offset=0, bad_cnt=0;
  
  char dis_name[256];


  if (argc!=3) {
    printf("Usage: %s <in_header_file> <out_cleaned_header_file>\n\n",argv[0]);
    printf("<in>\tName of input header file to clean");
    printf("<out>\tOutput name of cleaned header file");
    printf("\n\n");
    exit(1);
  }
  
  printf("\n\nDECODED SEASAT HEADER CLEANSING PROGRAM\n\n");
  printf("\topening files...\n");
  fpin=fopen(argv[1],"r");
  if (fpin==NULL) {printf("ERROR: Unable to open input file %s\n",argv[1]); exit(1);}
  fpout=fopen(argv[2],"w");
  if (fpout==NULL) {printf("ERROR: Unable to open output file %s\n",argv[2]); exit(1);}
  strcpy(dis_name,argv[2]);
  strcat(dis_name,".dis");
  fpdis=fopen(dis_name,"w");  /* destroy previous file if it existed */
  if (fpdis==NULL) {printf("ERROR: Unable to open output file %s\n",dis_name); exit(1);}
  fclose(fpdis);
  
  printf("\tinitializing histograms...\n");
  for (i=0;i<16;i++) {
    station_code_hist[i] = 0;   
    bits_per_sample_hist[i] = 0;
    prf_rate_code_hist[i] = 0;  
    lsd_year_hist[i] = 0;
  }
  for (i=0; i<MAX_DELAY; i++) delay_hist[i] = 0;
  for (i=0; i<400; i++) doy_hist[i] = 0;
  for (i=0; i<MAX_CLOCK_DRIFT; i++) clock_drift_hist[i] = 0;
  
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
    station_code_hist[hdr[i]->station_code]++;
    doy_hist[hdr[i]->day_of_year]++;
    clock_drift_hist[hdr[i]->clock_drift]++;
    delay_hist[hdr[i]->delay]++;
    lsd_year_hist[hdr[i]->lsd_year]++;
    bits_per_sample_hist[hdr[i]->bits_per_sample]++;
    prf_rate_code_hist[hdr[i]->prf_rate_code]++;
    
    times[i] = hdr[i]->msec;
    lines[i] = hdr[i]->major_cnt;
  }

  printf("\tgetting first median values\n");
  station_code_median    = get_median(station_code_hist,16);
  doy_median             = get_median(doy_hist,400);
  clock_drift_median     = get_median(clock_drift_hist,MAX_CLOCK_DRIFT);
  lsd_year_median        = get_median(lsd_year_hist,16);
  bits_per_sample_median = get_median(bits_per_sample_hist,16); 
  prf_rate_code_median   = get_median(prf_rate_code_hist,16);  
  delay_median           = get_median(delay_hist,MAX_DELAY); 

  printf("\t\tstation_code_median    = %i \n",station_code_median	);
  printf("\t\tdoy_median             = %li\n",doy_median		);
  printf("\t\tclock_drift_median     = %i \n",clock_drift_median	);
  printf("\t\tlsd_year_median        = %i \n",lsd_year_median	);
  printf("\t\tbits_per_sample_median = %i \n",bits_per_sample_median);
  printf("\t\tprf_rate_code_median   = %i \n",prf_rate_code_median  );
  printf("\t\tdelay_median           = %i \n",delay_median  	);

  // yax2bxc(lines,times,WINDOW_SIZE,&a,&b,&c);
  yaxb(lines,times,WINDOW_SIZE,&a,&b);
  
/*
  printf("============================================================================\
==================================================================\n"); 
  printf("BIN");
  printf("\t\t0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\n");
 
  printf("============================================================================\
==================================================================\n"); 
  printf("Station Code");
  for(i=0;i<16;i++) printf("\t%i",station_code_hist[i]);
  printf("\n");
      	    
  printf("============================================================================\
==================================================================\n"); 
  printf("Bits per Samp");
  for(i=0;i<16;i++) printf("\t%i",bits_per_sample_hist[i]);
  printf("\n");

  printf("============================================================================\
==================================================================\n"); 
  printf("PRF Rate Code");
  for(i=0;i<16;i++) printf("\t%i",prf_rate_code_hist[i]);
  printf("\n");
  
  printf("============================================================================\
==================================================================\n"); 
*/
/*
  printf("\nClock Drift Histogram values:\n");
  for(i=0;i<MAX_CLOCK_DRIFT;i++) if (clock_drift_hist[i]!=0) printf("%i\t%i\n",i,clock_drift_hist[i]);
*/


  /* Now, dump out the first WINDOW_SIZE/2 values to output file (corrected)
   ========================================================================*/
  printf("\tdumping initial lines to output file\n");
  for(i=0; i<WINDOW_SIZE/2; i++) {
  
    // double tmp = a*hdr[i]->major_cnt*hdr[i]->major_cnt+b*hdr[i]->major_cnt+c;
    double tmp = a*hdr[i]->major_cnt+b;
    if (fabs(tmp - hdr[i]->msec) > TOLERANCE) {
      printf("At %i bad value %li fixed value %lf\n",hdr[i]->major_cnt,hdr[i]->msec,tmp);
      msec = (long int) (tmp+0.5);
    } else msec = hdr[i]->msec;
	    
    fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    hdr[i]->major_cnt,
    hdr[i]->major_sync_loc,
    station_code_median,
    lsd_year_median,
    doy_median,
    msec,
    clock_drift_median,
    hdr[i]->no_scan_indicator_bit,
    bits_per_sample_median,
    hdr[i]->mfr_lock_bit,
    prf_rate_code_median,
    delay_median,
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
  
  /*=========================================================
    Loop through the rest of the file:
        - subtract one line from the histograms
	- read a new line and add it into the histograms,
	- get the histogram medians,
	- write a record out to file (corrected)
   =============================================================*/
  while (vals == 20) {
    if ((icnt%10000)==0) {printf("\tcleaning line %i\n",icnt);}
  
    /* Redo the time linear regression every so often 
     -----------------------------------------------*/
    if (icnt%(WINDOW_SIZE/20)==0) {
      // printf("Recalculating linear fit for MSEC now... ");
      for (i=0; i<WINDOW_SIZE; i++) {
        times[i] = hdr[i]->msec;
  	lines[i] = hdr[i]->major_cnt;
      }
      double old_a = a;
      double old_b = b;
      yaxb(lines,times,WINDOW_SIZE,&a,&b);
      if (a > 0.6073 || a < 0.607 ) { /* BAD do not use*/ a=old_a;b=old_b;/*printf("\tDiscarded\n");*/}
      else { /* printf("\tUsed\n"); */ offset = 0; }
    }
  
    /* remove the last values from the histograms 
     -------------------------------------------*/
    station_code_hist[hdr[curr]->station_code]--;
    doy_hist[hdr[curr]->day_of_year]--;
    clock_drift_hist[hdr[curr]->clock_drift]--;
    delay_hist[hdr[curr]->delay]--;
    lsd_year_hist[hdr[curr]->lsd_year]--;
    bits_per_sample_hist[hdr[curr]->bits_per_sample]--;
    prf_rate_code_hist[hdr[curr]->prf_rate_code]--;

    /* read and add in the next set of values 
     ---------------------------------------*/
    vals = get_values(fpin,hdr[curr]);
    if (vals == 20) {
      icnt++;
      station_code_hist[hdr[curr]->station_code]++;
      doy_hist[hdr[curr]->day_of_year]++;
      clock_drift_hist[hdr[curr]->clock_drift]++;
      delay_hist[hdr[curr]->delay]++;
      lsd_year_hist[hdr[curr]->lsd_year]++;
      bits_per_sample_hist[hdr[curr]->bits_per_sample]++;
      prf_rate_code_hist[hdr[curr]->prf_rate_code]++;

      /* get the new histogram median values 
       ------------------------------------*/
      station_code_median    = get_median(station_code_hist,16);
      doy_median             = get_median(doy_hist,400);
      clock_drift_median     = get_median(clock_drift_hist,MAX_CLOCK_DRIFT);
      lsd_year_median        = get_median(lsd_year_hist,16);
      bits_per_sample_median = get_median(bits_per_sample_hist,16); 
      prf_rate_code_median   = get_median(prf_rate_code_hist,16);  
      delay_median           = get_median(delay_hist,MAX_DELAY); 

      /* dump out the next line to the output file
       ------------------------------------------*/
      // double tmp = a*hdr[optr]->major_cnt*hdr[optr]->major_cnt+b*hdr[optr]->major_cnt+c;
      
      double tmp = a*(hdr[optr]->major_cnt+offset)+b;
      diff = fabs(tmp - hdr[optr]->msec);
      double sdiff = tmp - hdr[optr]->msec;
      
      if (diff > TOLERANCE) {
        printf("At %i bad value %li diff %lf cnt %i\n",hdr[optr]->major_cnt,hdr[optr]->msec,sdiff,bad_cnt);
        if (bad_cnt > SHIFT_GAP && (fabs(sdiff-save_diff)<0.9))
          {
           printf("\tSliding time window to fit possible discontinuity...last diff %lf this diff %lf\n",
      		save_diff,sdiff);
      	   int dir = -1.0*(tmp - hdr[optr]->msec)/diff;
      	   printf("\tdirection is %i\n",dir);
	   int save_offset = offset;
	   
	   if (diff > 4000.0) {
	     printf("ERROR: Unable to fix a gap of size %lf\n",diff);
	     printf("ERROR: This is equivalent to %lf seconds (%.0lf lines) of missing data!!!\n\n\n",
		     diff/1000.0,diff/1000.0*1647.0);
	     printf("Closing output files\n");
             fclose(fpout);
	     printf("Removing files %s and %s\n",argv[2],dis_name);
	     remove(argv[2]);
	     remove(dis_name);
	     exit(1);
	    }	   
	   while (diff > 1.0 && dir==1) {
	     offset+=dir;
	     tmp = a*(hdr[optr]->major_cnt+offset)+b;
	     diff = fabs(tmp - hdr[optr]->msec);
             printf("\tat offset %i: fixed value %lf diff %lf\n",offset,
	      hdr[optr]->major_cnt,hdr[optr]->msec,tmp,tmp-(double)hdr[optr]->msec);
           }
	   
	   if (dir==1 && offset > 5) {
	     printf("LOCATED DISCONTINUITY AT: %i ; LINES: %i \n",hdr[optr]->major_cnt,offset);
	     fpdis = fopen(dis_name,"a+");
	     fprintf(fpdis,"%i\t%i\t%i\t%lf\t%lf\n",hdr[optr]->major_cnt,save_offset,offset,a,b);
	     fclose(fpdis);
	   }
	   bad_cnt = 0; 
          } 
        else
          {
            if (bad_cnt==0) { save_diff = sdiff; bad_cnt++; }
            else { 
              if (fabs(sdiff-save_diff)<0.9) { save_diff=sdiff; bad_cnt++; } 
	      else { save_diff = sdiff; bad_cnt = 1; }
            }
          }
        msec = (long int) (tmp+0.5);
        // hdr[optr]->msec = msec; /* fix the past points for the next fit */
      } else { bad_cnt=0; msec = hdr[optr]->msec;}

      curr = (curr+1)%WINDOW_SIZE;
       
      fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
      hdr[optr]->major_cnt,
      hdr[optr]->major_sync_loc,
      station_code_median,
      lsd_year_median,
      doy_median,
      msec,
      clock_drift_median,
      hdr[optr]->no_scan_indicator_bit,
      bits_per_sample_median,
      hdr[optr]->mfr_lock_bit,
      prf_rate_code_median,
      delay_median,
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
  }

  /* Finally, dump out the last WINDOW_SIZE/2 values (corrected)
   =============================================================*/
  printf("\tdumping final lines to output file\n");
  for(i=0; i<WINDOW_SIZE/2; i++) {
    fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    hdr[optr]->major_cnt,
    hdr[optr]->major_sync_loc,
    station_code_median,
    lsd_year_median,
    doy_median,
    hdr[optr]->msec,
    clock_drift_median,
    hdr[optr]->no_scan_indicator_bit,
    bits_per_sample_median,
    hdr[optr]->mfr_lock_bit,
    prf_rate_code_median,
    delay_median,
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
  else printf("\n\nDone with calculations - read and wrote %i lines\n\n",icnt);
  
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

/*  printf("s->major_cnt);  %i\n",s->major_cnt);
    printf("s->major_sync_lo%li\n",s->major_sync_loc);
    printf("s->station_code)%i\n",s->station_code);
    printf("s->lsd_year);   %i\n",s->lsd_year);
    printf("s->day_of_year);%i\n",s->day_of_year);
    printf("s->msec);	    %li\n",s->msec);
    printf("s->clock_drift);%i\n",s->clock_drift);
    printf("s->no_scan_indic%i\n",s->no_scan_indicator_bit);
    printf("s->bits_per_samp%i\n",s->bits_per_sample);
    printf("s->mfr_lock_bit)%i\n",s->mfr_lock_bit);
    printf("s->prf_rate_code%i\n",s->prf_rate_code);
    printf("s->delay);	    %i\n",s->delay);
    printf("s->scu_bit);    %i\n",s->scu_bit);
    printf("s->sdf_bit);    %i\n",s->sdf_bit);
    printf("s->adc_bit);    %i\n",s->adc_bit);
    printf("s->time_gate_bit%i\n",s->time_gate_bit);
    printf("s->local_prf_bit%i\n",s->local_prf_bit);
    printf("s->auto_prf_bit)%i\n",s->auto_prf_bit);
    printf("s->prf_lock_bit)%i\n",s->prf_lock_bit);
    printf("s->local_delay_b%i\n",s->local_delay_bit);
*/    
    
  if(s->station_code    >=16) {
    // printf("ERROR: s->station_code    %i\n",s->station_code   ); 
    s->station_code=15;
  }
  if(s->day_of_year    >=400) {
    // printf("ERROR: s->day_of_year     %i\n",s->day_of_year    ); 
    s->day_of_year=399;
  }
  if(s->clock_drift>=MAX_CLOCK_DRIFT) {
    // printf("ERROR: s->clock_drift     %i\n",s->clock_drift    ); 
    s->clock_drift=MAX_CLOCK_DRIFT-1;
  }
  if(s->delay   >= MAX_DELAY) {
    // printf("ERROR: s->delay           %i\n",s->delay          ); 
    s->delay=MAX_DELAY-1;
  }
  if(s->lsd_year        >=16) {
    // printf("ERROR: s->lsd_year        %i\n",s->lsd_year       ); 
    s->lsd_year=15;
  }
  if(s->bits_per_sample >=16) {
    // printf("ERROR: s->bits_per_sample %i\n",s->bits_per_sample); 
    s->bits_per_sample=15;
  }
  if(s->prf_rate_code   >=16) {
    // printf("ERROR: s->prf_rate_code   %i\n",s->prf_rate_code  ); 
    s->prf_rate_code=15;
  }
    
  return(val);
}

int get_median(int *hist, int size) {
  int retval = -1, max = 0, i;
  for (i=0; i<size; i++) if (hist[i]>max) {max=hist[i]; retval=i;}
  if (retval==-1) { printf("Error getting histogram median value\n"); exit(1); }
  return(retval);    
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
   sum_x = 0.0;
   sum_xx = 0.0;
   sum_y = 0.0;
   sum_xy = 0.0;
   cnt = 0;
   
   for (i=0; i<n; i++)
    {
      if (x_vec[i] != 0.0) {
        sum_x += x_vec[i];
        sum_y += y_vec[i];
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

  // printf("Found coefficients y = %lf x + %lf\n",*a,*b);
   
}


/******************************************************************************
NAME:       yax2bxc.c

DESCRIPTION:    Computes a, b, and c for y = ax^2 + bx + c using quadratic
        	regression (least squares fit) given double vectors y and x
        	of length n.

PARAMETERS: 	x_vec   double[]     Input vector of X values
        	y_vec   double[]     Input vector of Y values
        	n   	int         Length of input vectors
        	a   	double*      Return x^2 coefficient
        	b   	double*      Return x coefficient
        	c   	double*      Return offset factor

HISTORY:        Ver 1.0   2/97   T. Logan    Initial Creation
                Ver 1.1   5/98   O. Lawlor   Changed to doubles
		Ver 1.2   9/12   T. Logan    Uses all doubles now
	`	
ALGORITHM REF:  Cheney, Ward & D. Kincaid, Numerical Mathematics and Computing,
            2nd Editn. pp 360-362. Brooks/Cole Pub. Co., Pacific Grove, Ca.
******************************************************************************/

void yax2bxc(double x_vec[],double y_vec[],int n,double *a,double *b,double *c)
{
 double x1, x2, x3, x4,     	/* Sum of x, x^2, x^3, x^4 */
       y1, yx, yx2;     	/* Sum of y, y*x, y*x^2    */
 double d1, d2, d3, d4, d5;     /* Intermediate Values     */
 double t1, t2, t3;     	/* Equation Solutions      */
 int   i;

 /* Calculate all of the first order sums */
 x1 = x2 = x3 = x4 = y1 = yx = yx2 = 0.0;

 for (i=0; i<n; i++) {
    x1  += x_vec[i];
    x2  += x_vec[i] * x_vec[i];
    x3  += x_vec[i] * x_vec[i] * x_vec[i];
    x4  += x_vec[i] * x_vec[i] * x_vec[i] * x_vec[i];
    y1  += y_vec[i];
    yx  += y_vec[i] * x_vec[i];
    yx2 += y_vec[i] * x_vec[i] * x_vec[i];
  }

 d1 = n*x2  - x1*x1;
 d2 = n*x3  - x1*x2;
 d3 = n*x4  - x2*x2;
 d4 = n*yx  - x1*y1;
 d5 = n*yx2 - x2*y1;

 t1 = (d1*d5 - d2*d4) / (d1*d3 - d2*d2);
 t2 = (d4 - d2*t1) / d1;
 t3 = (y1 - x2*t1 - x1*t2) / n;

 *a = t1;
 *b = t2;
 *c = t3;

 return;
}




