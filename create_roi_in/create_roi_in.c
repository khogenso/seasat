/******************************************************************************
NAME: create_roi_in - creates a ROI.in file from seasat HDR and state vectors
		      previously created by the ASF SEASAT PREP code.

SYNOPSIS: create_roi_in <infile>

DESCRIPTION:
	<infile> is a base name, assume that <infile>.dat and <infile>.hdr exist.

	- Read hdr file to get the start time and number of lines in the data segment
		- calculate the number of patches to process
	- Run dop.f (made into a subroutine) on the dat file
	- Fit output of doppler estimator with a 2nd order function
	- Find the correct state vector in the FIXED BODY state vectors
		- calculate spacecraft velocity
	- Run the state vector through get_peg_info (made into a subroutine)
	  to get the SCH Vel, SCH Acc, local earth radius, and spacecraft height
	- Create the <infile>.roi.in output file


EXTERNAL ASSOCIATES:
    NAME:               USAGE:
    ---------------------------------------------------------------

FILE REFERENCES:
    NAME:               USAGE:
    ---------------------------------------------------------------

PROGRAM HISTORY:
    VERS:   DATE:  AUTHOR:      PURPOSE:
    ---------------------------------------------------------------
    1.0	    4/12   T. Logan     Seasat Proof of Concept Project - ASF
    
HARDWARE/SOFTWARE LIMITATIONS:

ALGORITHM DESCRIPTION:

The order of outputs in the roi,in file is:

DESCRIPTION							VALUE
----------------------------------------------------------	-----------------------------------------------
First input data file					     <infile>.dat
Second input data file  				     /dev/null
Output data file					     <infile>.slc
Output amplitudes file  				     /dev/null
8lk output file 					     8lk
debug flag						     0
How many input bytes per line files 1 and 2		     13680 13680
How many good bytes per line, including header  	     13680 13680
First line to read  (start at 0)			     1
Enter # of range input patches  			     {Get nl from hdr values; patches = nl/11600}
First sample pair to use (start at zero)		     0
Azimuth Patch Size (Power of 2) 			     16384
Number of valid points in azimuth			     11600
Deskew the image					     n
Caltone % of sample rate				     0.25 0.25701904296875
Start range bin, number of range bins to process	     1 6840
Delta azimuth, range pixels for second file		     0 0
Image 1 Doppler centroid quad coefs (Hz/prf)		     {calculated from dop_est - 3 parameters}
Image 2 Doppler centroid quad coefs (Hz/prf)		     {copy above 3 values}
1 = use file 1 doppler, 2 = file 2, 3 = avg		     1
Earth Radius (m)					     {calculated from get_peg_info}
Body Fixed S/C velocities 1,2 (m/s)			     {calculate from EBEF state vector - 2 parameters}
Spacecraft height 1,2 (m)				     {calculated from get_peg_info - 2 parameters}
Planet GM						     0
Left, Right or Unknown Pointing 			     Right
SCH Velocity Vector 1					     {calculated from get_peg_info - 3 parameters}
SCH Velocity Vector 2					     {copy above 3 values}
SCH Acceleration Vector 1				     {calculated from get_pef_info - 3 parameters}
SCH Acceleration Vector 2				     {copy above 3 values}
Range of first sample in raw data file 1,2 (m)  	     {calculate from the hdr infomation - 2 values}
PRF 1,2 (pps)						     {calculate from the hdr information - 2 values}
i/q means, i1,q1, i2,q2 				     15.5 15.5 15.5 15.5
Flip i/q (y/n)  					     s
Desired azimuth resolution (m)  			     5  (what should this be???)
Number of azimuth looks 				     4  (what should this be???)
Range sampling rate (Hz)				     22765000
Chirp Slope (Hz/s)					     5.62130178e11
Pulse Duration (s)					     33.8e-6
Chirp extension points  				     0
Secondary range migration correction (y/n)		     n
Radar Wavelength (m)					     0.235
Range Spectral Weighting (1.=none, 0.54=Hamming)	     1.0
Fraction of range bandwidth to remove			     0 0
linear resampling coefs:  sloper, intr, slopea, inta	     0 0 0 0
linear resampling deltas: dsloper, dintr, dslopea, dinta     0 0 0 0
AGC file						     /dev/null
DWP file						     /dev/null or DWP file


ALGORITHM REFERENCES:

BUGS:

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "seasat.h"

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
int get_int_value(FILE *fp, const char token[], int *val, int from);
void estdop(FILE *fp, int sl, int nl, double *t1, double *t2, double *t3, double *iqmean);
void get_peg_info(double start_time, int nl, int prf, 
                  double *schvel, double *schacc, double *height, double *earthrad);
void spectra(FILE *fp,int sl, int nl,double iqmean,int *ocnt,double *ocal);

#define GOOD_SAMPLES  6840
#define GOOD_LINES    11600    
#define MAX_DWP_SHIFTS   20
#define DIGITIZATION_SHIFT 432
#define MAX_CALTONES  20

main(int argc, char *argv[])
{
  FILE *fpdat, *fphdr, *fproi, *fpstarthdr;
  char infile[256], outfile[256], hdrfile[256], starthdrfile[256], dwpfile[256];
  int err, nl, patches, prf;
  int from;
  double x,y,z,xdot,ydot,zdot,vel;
  double srf, time_length;
  double height, earthrad;
  double schvel[3], schacc[3];
  double t1, t2, t3;
  double iqmean;
  SEASAT_header_ext *hdr;
  SEASAT_header_ext *hdr1;

  double   start_sec,  current_sec, end_sec;
  double   time_from_start;
  int      start_date, current_date, end_date;
  int 	   start_year, current_year, end_year;
  int 	   end_sync;

  int tmp;
  double dtmp;
  double line_time_est;
  
  long int dwp_val[MAX_DWP_SHIFTS];
  long int dwp_line[MAX_DWP_SHIFTS];
  int      dwp_cnt, dwp_flag, dwp_min=65;
  
  int ncaltones;
  double caltones[MAX_CALTONES];

  julian_date s_date;
  julian_date e_date;
  ymd_date    s_ymd;
  ymd_date    e_ymd;
  hms_time    s_time;
  hms_time    e_time;
  
  int val, which;
  int i, start_line, end_line;

  if (argc!=2 && argc != 6) {
    printf("Usage: %s <infile_base_name> [-s <start_line> -e <end_line>]\n\n",argv[0]);
    printf("<infile_base_name>\tFile create roi input from. (assumes .dat and .hdr exist)\n");
    printf("\n\n");
    exit(1);
  }

  strcpy(infile,argv[1]); strcat(infile,".dat");
  strcpy(hdrfile,argv[1]); strcat(hdrfile,".hdr");

  start_line = 1;
  end_line = -99;
  
  if (argc == 6) {
    start_line = atoi(argv[3]);
    end_line = atoi(argv[5]);
  }

  if ((fphdr=fopen(hdrfile,"r"))==NULL) {printf("Error opening input file %s\n",hdrfile); exit(1);}
  hdr = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));
  hdr1 = (SEASAT_header_ext *) malloc(sizeof(SEASAT_header_ext));  
  
  printf("\n\n============================================================================\n");
  printf(" CREATING ROI.IN FILE FROM DATA %s\n",hdrfile);
  printf("============================================================================\n");

  for (i=0; i<MAX_DWP_SHIFTS; i++) {
    dwp_val[i] = 0;
    dwp_line[i] = 0;
  }

  /* hard coded the PRF to 1647 */
  prf = 1647;

  /* Create appropriate state vectors for this datatake
   ---------------------------------------------------*/
  val = get_values(fphdr, hdr);
  if (val!=20) {printf("ERROR: unable to read from header file\n"); exit(1);}
  
  s_date.year = 1970 + hdr->lsd_year;
  s_date.jd   = hdr->day_of_year;
  dtmp = (double) hdr->msec / 1000.0;
  date_sec2hms(dtmp,&s_time);
  
  /* set start year, day, second 
  -------------------------------------------------*/
  start_year = 1970 + hdr->lsd_year;
  start_date = hdr->day_of_year;
  start_sec  = (double) hdr->msec / 1000.0;
  
  /* seek to the start line the user requested
   ------------------------------------------*/
  for (i=1; i<start_line; i++) {
    val = get_values(fphdr, hdr);
    if (val!=20) {printf("ERROR: unable to read to specified start line in header file\n"); exit(1);}
  }
  current_year = 1970 + hdr->lsd_year;
  current_date = hdr->day_of_year;
  current_sec  = (double) hdr->msec / 1000.0;
  
  dwp_val[0]  = hdr->delay;
  dwp_line[0] = 0;
  dwp_cnt = 1;
  
  /* seek to the last line the user requested or else the end of file
   -------------------------------------------------------------------*/
  if (end_line == -99) {  /* read to the end of the file */
    which=0; nl = start_line-1;
    while (val==20) {
      nl++;
      if (which==0) { 
        val=get_values(fphdr,hdr1); 
	which=1; 
	if (dwp_val[dwp_cnt-1] != hdr1->delay) {
          dwp_val[dwp_cnt] = hdr1->delay;
          dwp_line[dwp_cnt] = nl-start_line;
          dwp_cnt++;
        }
      } else { 
        val=get_values(fphdr,hdr);  
	which=0; 
	if (dwp_val[dwp_cnt-1] != hdr->delay) {
          dwp_val[dwp_cnt] = hdr->delay;
          dwp_line[dwp_cnt] = nl-start_line;
          dwp_cnt++;
        }
      }
    }
    if (which==1) { /* we just read into hdr1, so hdr is good  */
      end_year = 1970+ hdr->lsd_year;
      end_date = hdr->day_of_year;
      end_sec  = (double) hdr->msec / 1000.0;
    } else {       /* we just read into hdr, so hdr1 is good */
      end_year = 1970+ hdr1->lsd_year;
      end_date = hdr1->day_of_year;
      end_sec  = (double) hdr1->msec / 1000.0;
    }
  } else {  /* just read to the end_line */
    nl = end_line-start_line+1;
    for (i=start_line; i<end_line; i++) {
      val = get_values(fphdr, hdr);
      if (val!=20) {printf("ERROR: unable to read to specified end line in header file\n"); exit(1);}
      if (dwp_val[dwp_cnt-1] != hdr->delay) {
        dwp_val[dwp_cnt] = hdr->delay;
        dwp_line[dwp_cnt] = i-start_line+1;
        dwp_cnt++;
      }
    }
    end_year = 1970+ hdr->lsd_year;
    end_date = hdr->day_of_year;
    end_sec  = (double) hdr->msec / 1000.0;
  }
   
  printf("Found start   time: %i %i %lf\n",start_year, start_date, start_sec);
  printf("Found current time: %i %i %lf\n",current_year, current_date, current_sec);
  printf("Found end     time: %i %i %lf\n",end_year, end_date, end_sec);
  printf("Found total lines : %i\n",nl);


  if (start_line != 1) {  /* we are processing a piece of the swath, add starting line to ROI.in file name */
    sprintf(dwpfile,"%s_line%i.dwp",argv[1],start_line);
  } else {
    strcpy(dwpfile,argv[1]); strcat(dwpfile,".dwp");
  }

  /* If we have at least one DWP change, need to create the DWP file */
  if (dwp_cnt > 1)  {  
    FILE *dwpfp = fopen(dwpfile,"w");
    int increasing;
    
    printf("Found DWP shifts in this scene; creating DWP file %s\n",dwpfile);
    for (i=0; i<dwp_cnt; i++) { if (dwp_val[i] < dwp_min) dwp_min = dwp_val[i]; }
    if (dwp_min == dwp_val[0]) increasing = 1;
    else increasing = 0;
    for (i=increasing; i<dwp_cnt; i++) {
        val = (dwp_val[i] - dwp_min) * DIGITIZATION_SHIFT;
	fprintf(dwpfp,"%i %i\n",dwp_line[i],val);
    }
    dwp_flag = 1;
  } else { dwp_flag = 0; dwp_min = dwp_val[0]; }

  printf("Propagating state vectors to requested time...\n");
  create_input_tle_file(s_date,s_time,"tle1.txt");
  propagate_state_vector("tle1.txt"); 
  printf("\n\nConverting state vectors from ECI to ECEF\n");
  fix_state_vectors(s_date.year,s_date.jd,s_time.hour,s_time.min,s_time.sec);
  
  // remove("tle1.txt");  
  // remove("propagated_state_vector.txt");

  /* Perform error checking on the times just read in 
  ---------------------------------------------------*/
  /* if years don't match, warning only - ignore it */
  // if (start_year != current_year) {printf("WARNING: Year of data take does not match!!!\n");}
    
  if (start_date != current_date) {
    if (current_date-start_date>1) { printf("ERROR: Bad current date found\n"); exit(1); }
    else current_sec += 86400.0;
  }
    
  if (end_date != current_date) {
    if (end_date-current_date>1) { printf("ERROR: Bad end date found\n"); exit(1); }
    else end_sec += 86400.0;
  }
   
  line_time_est = (double)nl / (double)prf;

  /* check for other timing errors */
  if (start_sec > end_sec) { 
      printf("WARNING: Data segment end time is before start of the datatake; fixing it (could be in error)!!!\n");
      end_sec = current_sec + line_time_est;
  }
    
  if (start_sec > current_sec) { printf("ERROR: Data segment time is before the start of the datatake!\n"); exit(1);}
  if (start_sec == 0) { printf("ERROR: Datatake time is ZERO!!!\n"); exit(1);}
  if (current_sec == 0) { printf("ERROR: Data segment start time is ZERO!!!\n"); exit(1);}
    
  time_from_start = current_sec - start_sec;
  time_length = end_sec - current_sec;
    
  if (fabs(line_time_est-time_length)>0.1) {
     printf("WARNING: Number of lines does not match time length\n"); 
     printf("WARNING: Time length from header: %lf; Line time estimate: %lf\n",time_length, line_time_est);
  }


/* Need to find the correct state vector for this data segment 
 ------------------------------------------------------------*/
  {
    FILE *fpvec;
    fpvec = fopen("fixed_state_vector.txt","r");
    int which = (int) (time_from_start+0.5);
    int i;
    double t;
    for (i=0; i<=which; i++)
      if (fscanf(fpvec,"%lf %lf %lf %lf %lf %lf %lf\n",&t,&x,&y,&z,&xdot,&ydot,&zdot)!=7) 
        { printf("ERROR: Unable to find state vector #%i in fixed_state_vector.txt file\n",which); exit(1); }
    fclose(fpvec);
    
    // remove("fixed_state_vector.txt");
  }

/* Calculate the slant range to the first pixel 
 ---------------------------------------------*/
  {
    double dwp, tau, pri;
    double c = 299792458.0;
    
    printf("found dwp min of %i\n",dwp_min);
    pri = 1.0 / (double)prf;
    dwp = ((double)dwp_min/64.0)*pri;
    tau = dwp + 9*pri;
    srf = tau * c / 2.0;
  }

/* Get the peg information that needed for ROI
 --------------------------------------------*/
  get_peg_info(time_from_start,nl,prf,schvel,schacc,&height, &earthrad);
  printf("Returned from get_peg_info\n");
 
/* Estimate the doppler centroid
 ------------------------------*/
  if ((fpdat=fopen(infile,"rb"))==NULL)  {printf("Error opening input file %s\n",infile); exit(1);}
  estdop(fpdat,start_line-1,nl,&t1,&t2,&t3,&iqmean);

/* Calculate the spectra and get the caltones
 -------------------------------------------*/
  fseek(fpdat,0,SEEK_SET);
  spectra(fpdat,start_line-1,nl,iqmean,&ncaltones,caltones);

/*=================================================================================
   NOW, ACTUALLY CREATE THE OUTPUT ROI FILE
 =================================================================================*/
  if (start_line != 1) {  /* we are processing a piece of the swath, add starting line to ROI.in file name */
    sprintf(outfile,"%s_line%i.roi.in",argv[1],start_line);
  } else {
    strcpy(outfile,argv[1]); 
    strcat(outfile,".roi.in");
  }
 
  if ((fproi=fopen(outfile,"w"))==NULL) {printf("Error opening output file %s\n",outfile); exit(1);}
  strcpy(outfile,argv[1]); strcat(outfile,".slc");

  printf("============================================================================\n");
  printf(" EMITTING FILE HEADER FILE NOW\n");
  printf("============================================================================\n");

/* First input data file */
  printf("First input data file: %s\n",infile);
  fprintf(fproi,"%s\n",infile);
  
/* Second input data file */
  printf("Second input data file: /dev/null\n");
  fprintf(fproi,"/dev/null\n");

/* Output data file */
  printf("Output data file: %s\n",outfile);
  fprintf(fproi,"%s\n",outfile);

/* Output amplitudes file */
  printf("Output amplitudes file: /dev/null\n");
  fprintf(fproi,"/dev/null\n");

/* 8lk output file */
  printf("8lk output file: 8lk\n");
  fprintf(fproi,"8lk\n");

/* debug flag */
  printf("debug flag: 0\n");
  fprintf(fproi,"0\n");
    
/* How many input bytes per line files 1 and 2 */
  printf("How many input bytes per line files 1 and 2: 13680 13680\n");
  fprintf(fproi,"13680 13680\n");

/* How many good bytes per line, including header */
  printf("How many good bytes per line, including header: 13680 13680\n");
  fprintf(fproi,"13680 13680\n");

/* First line to read  (start at 0) */
  printf("First line to read: %i\n",start_line);
  fprintf(fproi,"%i\n",start_line);

/* Enter # of range input patches */
  patches = nl / GOOD_LINES;
  if (patches ==0) patches = 1;
  printf("# of range input patches: %i\n",patches);
  fprintf(fproi,"%i\n",patches);
  
/* First sample pair to use (start at zero) */
  printf("First sample pair to use: 0\n");
  fprintf(fproi,"0\n");
  
/* Azimuth Patch Size (Power of 2) */
  printf("Azimuth Patch Size (Power of 2): 16384\n");
  fprintf(fproi,"16384\n");
  
/* Number of valid points in azimuth */
  printf("Number of valid points in azimuth: %i\n",GOOD_LINES);
  fprintf(fproi,"%i\n",GOOD_LINES);
  
/* Deskew the image */
  printf("Deskew the image: n\n");
  fprintf(fproi,"n\n");
  
/* Caltone % of sample rate */				
  printf("Number of caltones to remove: %i\n",ncaltones);
  fprintf(fproi,"%i\n",ncaltones);
  for (i=0;i<ncaltones;i++) {
    printf("\tCaltone %i: %lf\n",i,caltones[i]);
    fprintf(fproi,"%.14lf\n",caltones[i]);
  }
  
/* Start range bin, number of range bins to process */		
  printf("Start range bin, number of range bins to process: 1 %i\n",GOOD_SAMPLES);
  fprintf(fproi,"1 %i\n",GOOD_SAMPLES);
  
/* Delta azimuth, range pixels for second file */
  printf("Delta azimuth, range pixels for second file: 0 0\n");
  fprintf(fproi,"0 0\n");

/* Image 1 Doppler centroid quad coefs (Hz/prf)	*/
  printf("Doppler centroid quad coefs (Hz/prf): %lf %.8lf %.10lf\n",t1,t2,t3);
  fprintf(fproi,"%lf %.8lf %.12lf\n",t1,t2,t3);
  fprintf(fproi,"%lf %.8lf %.12lf\n",t1,t2,t3);

/* 1 = use file 1 doppler, 2 = file 2, 3 = avg */
  printf("1 = use file 1 doppler, 2 = file 2, 3 = avg: 1\n");
  fprintf(fproi,"1\n");
  
/* Earth Radius (m) */
  printf("Earth Radius of Curvature (m): %lf\n",earthrad);
  fprintf(fproi,"%lf\n",earthrad);

/* Body Fixed S/C velocities 1,2 (m/s) */
  vel = sqrt(xdot*xdot+ydot*ydot+zdot*zdot);
  printf("Body Fixed S/C velocities 1,2 (m/s): %lf %lf\n",vel,vel);
  fprintf(fproi,"%lf %lf\n",vel,vel);

/* Spacecraft height 1,2 (m) */
  printf("Spacecraft height 1,2 (m): %lf\n",height);
  fprintf(fproi,"%lf %lf\n",height,height);
  
/* Planet GM */
  printf("Planet GM: 0\n");
  fprintf(fproi,"0\n");
  
/* Left, Right or Unknown Pointing */
  printf("Left, Right or Unknown Pointing: Right\n");
  fprintf(fproi,"Right\n");
  
/* SCH Velocity Vector 1 & 2 */
  printf("SCH Velocity Vector 1 & 2: %lf %lf %lf\n",schvel[0],schvel[1],schvel[2]);
  fprintf(fproi,"%lf %lf %lf\n",schvel[0],schvel[1],schvel[2]);
  fprintf(fproi,"%lf %lf %lf\n",schvel[0],schvel[1],schvel[2]);
  
/* SCH Acceleration Vector 1 & 2 */
  printf("SCH Acceleration Vector 1 & 2: %lf %lf %lf\n",schacc[0],schacc[1],schacc[2]);
  fprintf(fproi,"%lf %lf %lf\n",schacc[0],schacc[1],schacc[2]);
  fprintf(fproi,"%lf %lf %lf\n",schacc[0],schacc[1],schacc[2]);
  
/* Range of first sample in raw data file 1,2 (m) */
  printf("Range of first sample in raw data file 1,2 (m): %lf %lf\n",srf,srf);
  fprintf(fproi,"%lf %lf\n",srf,srf);

/* PRF 1,2 (pps) */
  printf("PRF 1,2: %i %i\n",prf,prf);
  fprintf(fproi,"%i %i\n",prf,prf);
  
/* i/q means, i1,q1, i2,q2 */
  printf("i/q means, i1,q1, i2,q2: %lf %lf %lf %lf\n",iqmean,iqmean,iqmean,iqmean);
  fprintf(fproi,"%lf %lf %lf %lf\n",iqmean,iqmean,iqmean,iqmean);
  
/* Flip i/q (y/n) */
  printf("Flip i/q (y/n): s\n");
  fprintf(fproi,"s\n");
  
/* Desired azimuth resolution (m)  5  (what should this be???) */
  printf("Desired azimuth resolution (m): 5\n");
  fprintf(fproi,"5\n");

/* Number of azimuth looks         4  (what should this be???) */
  printf("Number of azimuth looks: 4\n");
  fprintf(fproi,"4\n");
  
/* Range sampling rate (Hz) */
  printf("Range sampling rate (Hz): 22765000\n");
  fprintf(fproi,"22765000\n");
  
/* Chirp Slope (Hz/s) */
  printf("Chirp Slope (Hz/s): 5.62130178e11\n");
  fprintf(fproi,"5.62130178e11\n");
  
/* Pulse Duration (s) */
  printf("Pulse Duration (s): 33.8e-6\n");
  fprintf(fproi,"33.8e-6\n");
  
/* Chirp extension points */
  printf("Chirp extension points: 0\n");
  fprintf(fproi,"0\n");
  
/* Secondary range migration correction (y/n) */
  printf("Secondary range migration correction (y/n): y\n");
  fprintf(fproi,"y\n");
  
/* Radar Wavelength (m) */					
  printf("Radar Wavelength (m): 0.235\n");
  fprintf(fproi,"0.235\n");
  
/* Range Spectral Weighting (1.=none, 0.54=Hamming) */
  printf("Range Spectral Weighting (1.=none, 0.54=Hamming): 1.0\n");
  fprintf(fproi,"1.0\n");
  
/* Fraction of range bandwidth to remove */
  printf("Fraction of range bandwidth to remove: 0 0\n");
  fprintf(fproi,"0 0\n");
  
/* linear resampling coefs:  sloper, intr, slopea, inta */
  printf("linear resampling coefs:  sloper, intr, slopea, inta: 0 0 0 0\n");
  fprintf(fproi,"0 0 0 0\n");
  
/* linear resampling deltas: dsloper, dintr, dslopea, dinta */
  printf("linear resampling deltas: dsloper, dintr, dslopea, dinta: 0 0 0 0\n");
  fprintf(fproi,"0 0 0 0\n");
  
/* AGC file */
  printf("AGC file: /dev/null\n");
  fprintf(fproi,"/dev/null\n");
  
/* DWP file */
  if (dwp_flag == 0) {
    printf("DWP file: /dev/null\n");
    fprintf(fproi,"/dev/null\n");
  } else {
    printf("DWP file: %s\n",dwpfile);
    fprintf(fproi,"%s\n",dwpfile);
  }

/*  fclose(fpdat); */
  fclose(fphdr);
  fclose(fproi);

  printf("============================================================================\n");
  printf(" CREATE_ROI_IN PROGRAM COMPLETED\n");
  printf("============================================================================\n\n\n");
}

int get_values(FILE *fp,SEASAT_header_ext *s)
{
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

