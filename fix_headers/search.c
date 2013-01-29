/******************************************************************************
NAME: discontinuity search - finds start of discontinuities, then applies
	fixes to header timings and inserts blank lines as needed to fill
	in discontinuous regions.

SYNOPSIS: fix_dis <indiscon> <indat> <inorghdr> <incleanhdr> <outdat> <outhdr>

DESCRIPTION:
	<indiscon> 	input discontinuity file
	<in>    	original input data and header files
	<incleanhdr> 	already cleaned header file to use values from
	<out> 		output data and header files with fill lines
	
This program follows the following algorithm:
    FOR each discontinuity from indiscon:
	read in RANGE lines before the discontinuity
	scan backwards from discontinuity to find actual start
	save result
    FOR each discontinuity start found:
    	read in/write out data and header until discontinuity start
	for length of discontinuity offset gap:
		write out header with corrected time and line
		insert blank line into output data file
	for length of run from start to previously found location:
		write out header with corrected time and line
		read in/write out data line
    WRITE out the rest of the file
    

EXTERNAL ASSOCIATES:
    NAME:               USAGE:
    ---------------------------------------------------------------

FILE REFERENCES:
    NAME:               USAGE:
    ---------------------------------------------------------------

PROGRAM HISTORY:
    VERS:   DATE:  AUTHOR:      PURPOSE:
    ---------------------------------------------------------------
    1.0	    11/12   T. Logan     Seasat Proof of Concept Project - ASF
    
HARDWARE/SOFTWARE LIMITATIONS:

ALGORITHM DESCRIPTION:

The order of columns in the header file is:

    major_cnt
    major_sync_loc,
    station_code,		
    lsd_year,			
    day_of_year,		
    msec,			
    clock_drift,		
    no_scan_indicator_bit,	
    bits_per_sample,		
    mfr_lock_bit,		
    prf_rate_code,		
    delay,			
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
void write_values(FILE *fp,SEASAT_header_ext *s);

#define RANGE 	     3000
#define SEARCH 	     200
#define MAX_DISCONS  1000
#define SAMPLES_PER_LINE 13680		/* decoded samples per output line */

main (int argc, char *argv[])
{
  double a[MAX_DISCONS];
  double b[MAX_DISCONS];
  int line[MAX_DISCONS];
  int start[MAX_DISCONS];
  int offset1[MAX_DISCONS];
  int offset2[MAX_DISCONS];
  
  int i, j, seek, val, dcnt = 0;
  int curr_line, cum_off;
  SEASAT_header_ext *hdr;
  long int tbuff[RANGE];
  long int lbuff[RANGE];
  FILE *fpdis, *fpin;
  FILE *fpin_dat, *fpin_hdr;
  FILE *fpout_dat, *fpout_hdr;
  int cnt, ptr;
  long int save;
  double tval;
  char cmd[256];
  unsigned char buf[SAMPLES_PER_LINE];
  int  total=0;
  
  char indat[256], inhdr[256];
  char outdat[256], outhdr[256], outdis[256];
  

  if (argc != 5) {
     printf("Usage: %s <indiscon> <in> <incleanhdr> <out>\n",argv[0]);
     printf("\tindiscon   - input discontinuity file\n");
     printf("\tin         - input data and header file base name\n");
     printf("\tincleanhdr - input cleaned header file\n");
     printf("\tout        - output data and header file base name\n\n");
     exit(1);
  }

  strcpy(indat,argv[2]); strcat(indat,".dat");
  strcpy(inhdr,argv[2]); strcat(inhdr,".hdr");
  strcpy(outdat,argv[4]); strcat(outdat,".dat");
  strcpy(outhdr,argv[4]); strcat(outhdr,".hdr");
  strcpy(outdis,argv[4]); strcat(outdis,".dis");

  hdr = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));

  fpdis = fopen(argv[1],"r");
  if (fpdis==NULL) { printf("WARNING: no discontinuity file found, assuming none needed\n"); }
  else { 
    while (fscanf(fpdis,"%i %i %i %lf %lf\n",&line[dcnt],&offset1[dcnt],&offset2[dcnt],&a[dcnt],&b[dcnt])==5) {
      printf("Read discontinuity %i\n",dcnt);
      printf("\tline   %i\n",line[dcnt]);
      printf("\tstart offset %i\n",offset1[dcnt]);
      printf("\tend offset %i\n",offset2[dcnt]);
      printf("\tCoeffs time = %lf line + %lf\n",a[dcnt],b[dcnt]);

      /* Read in RANGE lines before the discontinuity from the original (uncleaned) header file */  
      fpin = fopen(inhdr,"r");
      if (fpin == NULL) {printf("ERROR: Unable to open original input header file %s\n",inhdr); exit(1);}
      seek = line[dcnt] - RANGE;
      printf("\tseeking to line %i\n",seek);
      for (i=0;i<seek;i++) get_values(fpin,hdr);
      printf("\treading %i values\n",RANGE);
      for (i=0;i<RANGE;i++) {
        get_values(fpin,hdr);
	lbuff[i] = hdr->major_cnt;
	tbuff[i] = hdr->msec;
      }
      
      /* Scan backwards to find the start of this discontinuity */
      cnt = 0;
      ptr = RANGE;
      while (cnt < SEARCH && ptr > 0) {
          ptr--;
	  tval = a[dcnt]*(lbuff[ptr]+offset2[dcnt])+b[dcnt];
	  // printf("At %i, time=%li calc=%lf diff %lf cnt %i\n",
	  //  seek+ptr,tbuff[ptr],tval,fabs(tval-tbuff[ptr]),cnt);
	  if (fabs(tval-tbuff[ptr])<1.5) {
	    cnt = 0;
	    save = seek+ptr+1;
	  } else { cnt++; }
      }

      if (ptr <= 0) { 
        printf("ERROR: UNABLE TO FIX THIS DISCONTINUITY - it is too long\n");
	printf("ERROR: Probably need to increase the RANGE value in this code\n");
	exit(1);
      }
      
      printf("\tDISCONTINUITY #%i: start is %i\n",dcnt,save);
      start[dcnt]=save;
      fclose(fpin);
      dcnt++;
    }
    fclose(fpdis);
  }
  printf("Finished finding discontinuity starting points...\n");
  printf("Applying results to data and header files...\n");
  
  /* Now apply each of the discontinuities to the .dat and .hdr files */
  fpin_hdr = fopen(argv[3],"r");
  if (fpin_hdr == NULL) {printf("ERROR: Unable to open input cleaned header file %s\n",argv[3]); exit(1);}
  fpin_dat = fopen(indat,"rb");
  if (fpin_dat == NULL) {printf("ERROR: Unable to open input data file %s\n",indat); exit(1);}  
  fpout_hdr = fopen(outhdr,"w");
  if (fpout_hdr == NULL) {printf("ERROR: Unable to open output header file %s\n",outhdr); exit(1);}
  fpout_dat = fopen(outdat,"wb");
  if (fpout_dat == NULL) {printf("ERROR: Unable to open output data file %s\n",outdat); exit(1);}  
  
  if (dcnt>0) {
    fpdis = fopen(outdis,"w");
    if (fpdis==NULL) { printf("ERROR: Unable to open output dis file %s\n",outdis); exit(1);}
    fprintf(fpdis,"LINE\tGAP\n");
  }

  curr_line = 0; cum_off = 0;
  
  /* Deal with the discontinuities one at a time */
  for (i=0; i<dcnt; i++) {
    /* read in and write out lines until discontinuity is hit */
    for (j=curr_line; j<start[i]+cum_off; j++) {
      fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
      fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); total++;
      get_values(fpin_hdr,hdr);
      hdr->major_cnt=j;
      write_values(fpout_hdr,hdr);
    }
    curr_line = j;
    printf("\twrote unchanged to line %i\n",curr_line);
    for (j=0; j<SAMPLES_PER_LINE; j++) buf[j] = 0;
    
    /* repeat the header with correct times and insert blanks for length of gap */
    printf("\tfilling in a gap of %i\n",offset2[i]-offset1[i]);
    for (j=curr_line; j<curr_line+(offset2[i]-offset1[i]); j++) {
      fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); total++;
      hdr->msec = a[i]*(j-cum_off+offset1[i]) + b[i];  /* correct for fact that a,b are referenced to original lines */
      hdr->major_cnt = j;
      write_values(fpout_hdr,hdr);
    }
    curr_line = j;
    printf("\twrote fill values to line %i\n",curr_line);
    printf("\tchanging times for next %i lines\n",line[i]-start[i]);
   
    /* read in and write out the rest of this discontinuity fixing lines and times as we go */   
    for (j=curr_line; j<curr_line+(line[i]-start[i]); j++)  {
      fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
      fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); total++;
      get_values(fpin_hdr,hdr);
      hdr->major_cnt=j;
      hdr->msec = a[i]*(j-cum_off+offset1[i]) + b[i];
      write_values(fpout_hdr,hdr);
    }
    curr_line = j;
    fprintf(fpdis,"%i\t%i\n",start[i]+cum_off,offset2[i]-offset1[i]);
    cum_off = cum_off + offset2[i] - offset1[i];
    printf("\twrote fixed values to line %i (formerly %i)\n",curr_line,curr_line-cum_off);
    
  }
  
  if (dcnt>0) {fclose(fpdis);}
  
  /* Deal with the rest of the file */
  printf("Done with discontinuities, reading/writing rest of the file\n");  
  val = get_values(fpin_hdr,hdr);
  while (val==20) {
    fread(buf,SAMPLES_PER_LINE,1,fpin_dat);
    fwrite(buf,SAMPLES_PER_LINE,1,fpout_dat); total++;
    hdr->major_cnt = curr_line;
    write_values(fpout_hdr,hdr);
    val = get_values(fpin_hdr,hdr);
    curr_line++;
  }
  
  printf("\n");
  printf("Done correcting file, wrote %i lines of output\n\n",total);

  fclose(fpin_dat);
  fclose(fpin_hdr);
  fclose(fpout_dat);
  fclose(fpout_hdr);

  exit(0);
}


int get_values(FILE *fp,SEASAT_header_ext *s) {
  int val;
  if (s==NULL) {printf("empty pointer passed to get_values\n"); exit(1);}
  val = fscanf(fp,"%i %li %i %i %i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    &(s->major_cnt),&(s->major_sync_loc),&(s->station_code),&(s->lsd_year),
    &(s->day_of_year),&(s->msec),&(s->clock_drift),&(s->no_scan_indicator_bit),
    &(s->bits_per_sample),&(s->mfr_lock_bit),&(s->prf_rate_code),&(s->delay),
    &(s->scu_bit),&(s->sdf_bit),&(s->adc_bit),&(s->time_gate_bit),&(s->local_prf_bit),
    &(s->auto_prf_bit),&(s->prf_lock_bit),&(s->local_delay_bit));
  return(val);
}

void write_values(FILE *fp,SEASAT_header_ext *s) {
  if (s==NULL) {printf("empty pointer passed to write values\n"); exit(1);}
  if (fp==NULL) {printf("null file pointer passed to write values\n"); exit(1);}
  fprintf(fp,"%i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    s->major_cnt,s->major_sync_loc,s->station_code,s->lsd_year,
    s->day_of_year,s->msec,s->clock_drift,s->no_scan_indicator_bit,
    s->bits_per_sample,s->mfr_lock_bit,s->prf_rate_code,s->delay,
    s->scu_bit,s->sdf_bit,s->adc_bit,s->time_gate_bit,s->local_prf_bit,
    s->auto_prf_bit,s->prf_lock_bit,s->local_delay_bit);
}

	
