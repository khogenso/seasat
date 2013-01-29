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
#define TOLERANCE  257		/* how far off a time value can be from estimate */


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
void copy_hdr(SEASAT_header_ext *in,SEASAT_header_ext *out);

main(int argc, char *argv[])
{
  FILE *fpin,*fpout;
  SEASAT_header_ext **hdr;
  int i, j, k, vals;
  int icnt = 0, ocnt = 0, curr = 0, optr = WINDOW_SIZE/2;
  int fcnt = 0, n=0;
 
  long int msec;
  double times[WINDOW_SIZE];
  double lines[WINDOW_SIZE];
  double a, b, c;
  double diff, save_diff, tmp;
  double pri;
  
  int offset=0, bad_cnt=0;

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
  
  /* Allocate the array of headers */
  hdr = (SEASAT_header_ext **) malloc(sizeof(SEASAT_header_ext *)*WINDOW_SIZE);
  if (hdr == NULL) {printf("ERROR: unable to allocate hdr\n"); exit(1);}
  for (i=0; i<WINDOW_SIZE; i++)
    hdr[i] = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));
  for (i=0; i<WINDOW_SIZE; i++)
    if (hdr[i]==NULL) {printf("Error in malloc\n"); exit(1);}

  vals = get_values(fpin,hdr[1]);
  if (vals != 20) {printf("ERROR: can't read from input file\n"); exit(1);}
  icnt++;

  pri = 1000* (1.0 / 1647.0);

  while (vals == 20) {
    if ((icnt%10000)==0) {printf("\tcleaning line %i\n",icnt);}

    copy_hdr(hdr[1],hdr[0]);
    fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
      hdr[1]->major_cnt,hdr[1]->major_sync_loc,hdr[1]->station_code,
      hdr[1]->lsd_year,hdr[1]->day_of_year,hdr[1]->msec,hdr[1]->clock_drift,
      hdr[1]->no_scan_indicator_bit,hdr[1]->bits_per_sample,hdr[1]->mfr_lock_bit,
      hdr[1]->prf_rate_code,hdr[1]->delay,hdr[1]->scu_bit,hdr[1]->sdf_bit,
      hdr[1]->adc_bit,hdr[1]->time_gate_bit,hdr[1]->local_prf_bit,hdr[1]->auto_prf_bit,
      hdr[1]->prf_lock_bit,hdr[1]->local_delay_bit);
    ocnt++;

    n = 1;
    vals = get_values(fpin,hdr[1]);
    if (vals==20) icnt++;
    
    while (vals==20 && hdr[1]->msec == hdr[0]->msec) {
       double offset = pri*(double)n;
       hdr[1]->msec = hdr[0]->msec + (int) offset;
       // printf("Fixing... was %li now %li offset %lf\n",hdr[0]->msec,hdr[1]->msec,offset);
       fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
         hdr[1]->major_cnt,hdr[1]->major_sync_loc,hdr[1]->station_code,
         hdr[1]->lsd_year,hdr[1]->day_of_year,hdr[1]->msec,hdr[1]->clock_drift,
         hdr[1]->no_scan_indicator_bit,hdr[1]->bits_per_sample,hdr[1]->mfr_lock_bit,
         hdr[1]->prf_rate_code,hdr[1]->delay,hdr[1]->scu_bit,hdr[1]->sdf_bit,
         hdr[1]->adc_bit,hdr[1]->time_gate_bit,hdr[1]->local_prf_bit,hdr[1]->auto_prf_bit,
         hdr[1]->prf_lock_bit,hdr[1]->local_delay_bit);
       ocnt++;
       n++;
       vals = get_values(fpin,hdr[1]);
       if (vals==20) icnt++;
       if ((int)offset>=1) fcnt++;
    }
  }
       
 /*
  fprintf(fpout,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    hdr[1]->major_cnt,hdr[1]->major_sync_loc,hdr[1]->station_code,
    hdr[1]->lsd_year,hdr[1]->day_of_year,hdr[1]->msec,hdr[1]->clock_drift,
    hdr[1]->no_scan_indicator_bit,hdr[1]->bits_per_sample,hdr[1]->mfr_lock_bit,
    hdr[1]->prf_rate_code,hdr[1]->delay,hdr[1]->scu_bit,hdr[1]->sdf_bit,
    hdr[1]->adc_bit,hdr[1]->time_gate_bit,hdr[1]->local_prf_bit,hdr[1]->auto_prf_bit,
    hdr[1]->prf_lock_bit,hdr[1]->local_delay_bit);
  ocnt++;
  */
       
  if (icnt != ocnt) printf("ERROR: input/output don't match; read %i wrote %i\n",icnt,ocnt);
  else printf("\n\nDone with calculations - read and wrote %i lines; fixed %i time values (%f%%)\n\n",
  	icnt,fcnt, 100.0*(float)fcnt/(float)icnt);
  
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

void copy_hdr(SEASAT_header_ext *in,SEASAT_header_ext *out)
{
    out->major_cnt = 		  in->major_cnt ; 
    out->major_sync_loc = 	  in->major_sync_loc ; 
    out->station_code = 	  in->station_code ; 
    out->lsd_year = 		  in->lsd_year ; 
    out->day_of_year =  	  in->day_of_year ; 
    out->msec = 		  in->msec ; 
    out->clock_drift =  	  in->clock_drift ; 
    out->no_scan_indicator_bit =  in->no_scan_indicator_bit ; 
    out->bits_per_sample = 	  in->bits_per_sample ; 
    out->mfr_lock_bit = 	  in->mfr_lock_bit ; 
    out->prf_rate_code = 	  in->prf_rate_code ; 
    out->delay = 		  in->delay ; 
    out->scu_bit = 		  in->scu_bit ; 
    out->sdf_bit = 		  in->sdf_bit ; 
    out->adc_bit = 		  in->adc_bit ; 
    out->time_gate_bit = 	  in->time_gate_bit ; 
    out->local_prf_bit = 	  in->local_prf_bit ; 
    out->auto_prf_bit = 	  in->auto_prf_bit ; 
    out->prf_lock_bit = 	  in->prf_lock_bit ; 
    out->local_delay_bit = 	  in->local_delay_bit ; 
}
