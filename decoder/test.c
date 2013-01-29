/******************************************************************************
NAME: SeasatPrep - preps seasat format raw data into L0 raw files to be 
		   processed by ROI or similar SAR processor

SYNOPSIS:  SeasatPrep <infile> <outfile>

DESCRIPTION:

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

ALGORITHM REFERENCES:

BUGS:

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seasat.h"

SEASAT_raw_header r;
SEASAT_header s;
SEASAT_aux a;

#define MAX_FILE_LENGTH 97600   /* maximum number of lines to write to a data segment */
				/* estimated at 8 patches * 11,600 + (16k-11,600)     */

main(int argc,char *argv[])
{
  unsigned int testval, newval;
  unsigned char *buf;
  unsigned char inbuf[INT_FRAME_LEN+1];
  char inname[256], outname[256], outheadername[256];
  int  i, optr=0;
  int  headers=1000;
  FILE *fpin, *fpout=NULL, *fp_hdr=NULL;
  FILE *frame_file1, *frame_file2;
  FILE *fptmp, *fp_time, *fp_all_hdr;

  long seeking = 0;			/* keep track of how far we seek to find a sync */
  int found_cnt = 0;			/* number of sync found */
  int major_cnt = 0;			/* number of major frames found */
  int this_major_cnt = 0;		/* number of major frames found in this dataset */
  int pre_cnt = 0;			/* number of syncs before first dataset found */
  int fill_cnt = 0;			/* number of fill data flags found */
  int good_cnt = 0;			/* number of non-fill data flags found */
  int done = 0;				/* set when we are within one frame of the end of file */
  int lock = 0;				/* set when we have a lock on a dataset */
  int max_consecutive_fills=0;		/* how many fills did we find in a row at most? */
  long file_length;
  long offset;
  long last_sync=0;
  long this_sync=0;
  long major_sync_loc;			/* file location of last major sync */
  int  major_sync_num;			/* sync number of last major sync */
  unsigned char rtmp1,rtmp2, tmp, tmpc;
  unsigned char obuff[SAMPLES_PER_LINE+SAMPLES_PER_FRAME*10];
  double dtmp;
  int this_frame, next_frame, last_frame=-1;
  int expected_frame_no=-1;
  int frames_this_line=-1;
  int last_max_frame_no=-1;
  int fixed_frame_nos=0;
  int non_fixed_frame_nos=0;
  int l2_frame=-1, l3_frame=-1;		 /* last last and last last last frame numbers     */
  int dataset_number=0;			 /* which dataset we are working on, starts at 000 */
  long int last_bad=0;			 /* sync location for last sync that was fill data */
  int contiguous_miss=0;		 /* number of contiguous fill data frames read     */
  int partial_line_cnt=0;		 /* how many decoded lines are partials            */
  unsigned char end_of_dataset=0;	 /* set when number of contiguous misses is greater than MAX_CONTIGUOUS_MISSES
  					    signifies the end of one file and the (potential) start of another file */
  unsigned char DUMP_FIXED_FRAMES=0;	 /* set to 1 - all of the fixed frame numbers dumped to file	*/
  unsigned char DUMP_NON_FIXED_FRAMES=0; /* set to 1 - all of the original frame numbers dumped to file */
  unsigned char DISPLAY_FRAME_FIXES=0;   /* set to 1 to see which frames numbers have been fixed	*/
  unsigned char DUMP_TIMES=0;		 /* set to 1 to see all of the decoded times 			*/
  unsigned char DUMP_ALL_HEADERS=0;      /* set to 1 to see all of the decoded header values 		*/
  unsigned char DUMP_ALL_SYNCS=0;        /* set to 1 to see all of the sync locations 			*/
  int error = 0;

  int read_cnt;
  int aligned;

  julian_date start_date;
 
  if (argc != 3) {
    printf("Usage: %s <inname> <outname>\n",argv[0]);
    printf("\n");
    printf("inname  \tname of input RAW Seasat file (sync'd)\n");
    printf("outname \tbase name of output RAW Seasat file (prep'd)\n\n");
    exit(1);
  }
  strcpy(inname,argv[1]);

  r.aux0 =(SEASAT_aux_0 *) malloc(sizeof(SEASAT_aux_0));
  r.aux1 =(SEASAT_aux_1 *) malloc(sizeof(SEASAT_aux_1));
  r.aux2 =(SEASAT_aux_2 *) malloc(sizeof(SEASAT_aux_2));
  r.aux3 =(SEASAT_aux_3 *) malloc(sizeof(SEASAT_aux_3));
  r.aux4 =(SEASAT_aux_4 *) malloc(sizeof(SEASAT_aux_4));
  r.aux5 =(SEASAT_aux_5 *) malloc(sizeof(SEASAT_aux_5));
  r.aux6 =(SEASAT_aux_6 *) malloc(sizeof(SEASAT_aux_6));
  r.aux7 =(SEASAT_aux_7 *) malloc(sizeof(SEASAT_aux_7));
  r.aux8 =(SEASAT_aux_8 *) malloc(sizeof(SEASAT_aux_8));
  r.aux9 =(SEASAT_aux_9 *) malloc(sizeof(SEASAT_aux_9));

  for (i=0; i<SAMPLES_PER_LINE; i++) obuff[i] = 0;
  
  /* open input file and determine file length */
  fpin = fopen(inname,"rb");
  if (fpin==NULL) { printf("ERROR: unable to open input file %s\n",inname); exit(1); }
  fseek(fpin,0L,SEEK_END);
  file_length = ftell(fpin);
  fseek(fpin,0L,SEEK_SET);
  printf("Total file length is %li\n",file_length);
  printf("Starting at %li\n",ftell(fpin));
  
  if (DUMP_FIXED_FRAMES==1) {
    frame_file1 = fopen("frames_fixed.out","w");
    if (frame_file1==NULL) { printf("ERROR: unable to open output file frames_fixed.out\n"); exit(1);}
  }
  if (DUMP_NON_FIXED_FRAMES==1) {
    frame_file2 = fopen("frames_non_fixed.out","w");
    if (frame_file2==NULL) { printf("ERROR: unable to open output file frames_non_fixed.out\n"); exit(1);}
  }
  if (DUMP_TIMES==1) {
    fp_time = fopen("times.out","w");
    if (fp_time==NULL) { printf("ERROR: unable to open output file times.out\n"); exit(1);}
  }
  if (DUMP_ALL_HEADERS==1) {
    sprintf(outname,"%s.headers",argv[2]);
    fp_all_hdr = fopen(outname,"w");
    if (fp_all_hdr==NULL) { printf("ERROR: unable to open output file %s\n"); exit(1);}
  }
  if (DUMP_ALL_SYNCS==1) {  
    fptmp=fopen("sync_loc.txt","w");
    if (fptmp==NULL) { printf("ERROR: unable to open output file sync_loc.txt\n"); exit(1);}
  }

  read_cnt = 0;
  aligned = -1;

  /*===========================================================================================
                                           START MAIN LOOP 
  ===========================================================================================*/
  while (!feof(fpin)&&!done&&!error) {
    long where;
    this_sync = ftell(fpin);
    
    if (read_cnt == 28434) 		/* special case - only 147 byte frame here */
      { aligned = 2; read_cnt = 0;}
    else {
      if (read_cnt%2 == 0) 		/* read in 147.5 byte aligned frame */
        { aligned = 1; read_cnt++; }
      else 		 		/* read in 147.5 byte non-aligned frame */
        { aligned = 0; read_cnt++; }
    }

    fread(inbuf,INT_FRAME_LEN*sizeof(unsigned char),1,fpin);
    if (aligned==1 || aligned==2) fseek(fpin,-1L,SEEK_CUR);

    found_cnt++;

//    if (found_cnt>1 && this_frame != 127) {
    if (found_cnt>1) {
        l3_frame = l2_frame;
	l2_frame = last_frame;
      	last_frame = this_frame;
    }

    /* If we have an aligned sync code, do this to get the frame number */
    if (aligned==1 || aligned==2) {
        tmpc = inbuf[3];
	a.fill_flag = tmpc >> 7;
	a.frame_no = tmpc & 0177;
    }

    /* If we have an unaligned sync code, do this to get the frame number */  	
    if (aligned==0) {
        rtmp1 = inbuf[3];
	rtmp2 = inbuf[4];
	tmpc = rtmp1 << 4 | rtmp2 >> 4;
	a.fill_flag = tmpc >> 7;
	a.frame_no = tmpc & 0177;
    }
      
    this_frame = a.frame_no;
    next_frame = get_next_frameno(fpin,aligned);

    /* Check for bit errors in frame no 0 as this causes BIG problems 
     ---------------------------------------------------------------*/
    if ((this_frame-last_frame)!=1 && this_frame != 0 && (a.frame_no!=127 || next_frame!=127)) {  

	/* if the next frame is 1 and last was either 59 or 60, 
	   assume this is frame zero... might be a bad assumption???  
	  ---------------------------------------------------------*/
	if (next_frame==1 && (last_frame==59 || last_frame==60)) { 
	    if (DISPLAY_FRAME_FIXES==1)
	        printf("1: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO 000\n",
	    	l3_frame,l2_frame,last_frame,this_frame,next_frame);
	    this_frame = 0; fixed_frame_nos++;
	}
	
	/* if  (next_frame-last_frame)==2, put this_frame in sequence 
	 -----------------------------------------------------------*/
	else if (next_frame - last_frame == 2) {
	    if (DISPLAY_FRAME_FIXES==1)
	        printf("2: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO %3.3i\n",
	    	l3_frame,l2_frame,last_frame,this_frame,next_frame,last_frame+1);
	  this_frame = last_frame+1; fixed_frame_nos++;
	}
	
	/* if last two frames are in sequence, put this one in sequence
	 -------------------------------------------------------------*/
	else if (((l2_frame-l3_frame==1) || (l2_frame-l3_frame)%59 == 0 || (l2_frame-l3_frame)%60 == 0)
		 && (last_frame-l2_frame==1)) {
	  /* - can't ever have two lines in a row that have 60 minor frames
	     - otherwise, we are good to set this frame into sequence */
	  if (last_max_frame_no<60 || last_frame < 59) { 
    	    if (DISPLAY_FRAME_FIXES==1)
	        printf("3: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO %3.3i\n",
	    	l3_frame,l2_frame,last_frame,this_frame,next_frame,last_frame+1);
	    this_frame=last_frame+1; fixed_frame_nos++;
	  } else { 
	    /* If the last frame was either 59 or 60, then this HAS to be frame zero
	     ------------------------------------------------------------------------*/
    	    if (DISPLAY_FRAME_FIXES==1)
	        printf("4: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO 000\n",
	    	l3_frame,l2_frame,last_frame,this_frame,next_frame);
	    this_frame=0; fixed_frame_nos++;
	  }
	} 
	
	/* if last2 and last3 frames are in sequence and last frame 
	   is ZERO, then set this frame to 1
	 ---------------------------------------------------------*/
	else if ((l2_frame-l3_frame==1) && (last_frame==0)) {
	  if (DISPLAY_FRAME_FIXES==1)
            printf("5: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO 001\n",
	    	l3_frame,l2_frame,last_frame,this_frame,next_frame);
	  this_frame = 1; fixed_frame_nos++;
	}
	
	/* If we got to here, then we did not fix the error
	 -------------------------------------------------*/
	else if (lock==1) {
	    printf("6: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; NOT FIXED!! (Line %i)\n",
	    	l3_frame,l2_frame,last_frame,this_frame,next_frame,major_cnt);
	    non_fixed_frame_nos++;
	}
    }
     
    /* Check for bad frame zeros!  This can only be frame ZERO if the
       last frame was 59 or greater and the next frame is 1.
       -- Update - this condition doesn't allow 1/2 lines to exist,
          but they do exist in the raw files!!!
     -------------------------------------------------------------- */
    if (lock==1 && this_frame==0 && abs(this_frame-last_frame)<59 && next_frame!=1) {
      if (DISPLAY_FRAME_FIXES==1)
   	  printf("7: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO %3.3i\n",
   	      l3_frame,l2_frame,last_frame,this_frame,next_frame,last_frame+1);
      this_frame = last_frame+1; fixed_frame_nos++;
    }

    if (lock==1) frames_this_line++;

    /* We can NEVER have a frame number greater than 60!!!
     ----------------------------------------------------*/
    if (lock==1 && this_frame > 60 && frames_this_line>59) {
      if (DISPLAY_FRAME_FIXES==1)
	    printf("8: L3=%3.3i L2=%3.3i Last=%3.3i This=%3.3i Next=%3.3i; FIXED TO 000\n",
	    l3_frame,l2_frame,last_frame,this_frame,next_frame);
	this_frame=0; fixed_frame_nos++;
    }        

    /* printf("frame no %i; fill flag %i\n",a.frame_no, a.fill_flag); */
    if (DUMP_FIXED_FRAMES==1) { if (this_frame!=0 && lock==1) fprintf(frame_file1,"%2.2i ",this_frame); }
    if (DUMP_NON_FIXED_FRAMES==1) { if (this_frame!=0) fprintf(frame_file2,"%2.2i ",a.frame_no); }
  
    /* set pointer to data portion of the frame, shift if necessary */
    if (aligned==1 || aligned==2) buf = &(inbuf[4]); 
    else /* if (aligned==0) */ { 
      buf = &(inbuf[5]);
      shift_buffer(buf,rtmp2); 
    }
      
    /* start of a new major frame */
    if (this_frame == 0) { 
	headers = 0; 
	major_cnt++;
	this_major_cnt++;
	long int last_major_sync_loc = major_sync_loc;
	major_sync_loc = ftell(fpin)-(INT_FRAME_LEN-1); 
	if (aligned==0) major_sync_loc--;
	major_sync_num = found_cnt;
	if (lock == 0) {
          printf("==========================================================================\n");
  	  printf(" S Y N C    L O C K    E S T A B L I S H E D at sync %li, byte %li\n",major_sync_num,major_sync_loc);
          printf("==========================================================================\n");
          lock = 1;
          end_of_dataset=0;
	  last_frame = -1;
	  l2_frame = -1;
	  l3_frame = -1;
	} else {
  	  last_max_frame_no = last_frame;
  	  if (last_max_frame_no != 59 && last_max_frame_no != 60) {
	    printf("ERROR: Found range line #%li with %i frames in it!!!\n",this_major_cnt-1,last_max_frame_no);
	    partial_line_cnt++;
	  }
	}
	
	if (DUMP_FIXED_FRAMES==1) {fprintf(frame_file1,"\n"); fprintf(frame_file1,"%2.2i ",this_frame);}
	if (DUMP_NON_FIXED_FRAMES==1) {fprintf(frame_file2,"\n"); fprintf(frame_file2,"%2.2i ",a.frame_no);}
	frames_this_line = 1;
	
	if (this_major_cnt == 1) {
	  sprintf(outname,"%s_%3.3i.dat",argv[2],dataset_number);
      	  sprintf(outheadername,"%s_%3.3i.hdr",argv[2],dataset_number);
	  printf("Opening file %s for output\n",outname);
	  if (fpout==NULL) fpout = fopen(outname,"wb");
	  if (fpout==NULL) { printf("ERROR: unable to open output file %s\n",outname); exit(1); }
	  if (fp_hdr==NULL) fp_hdr = fopen(outheadername,"w");
	  if (fp_hdr==NULL) { printf("ERROR: unable to open output file %s\n",outheadername); exit(1); }
	  dataset_number++;
        }

	/* DUMP OUT LAST BUFFER OF DECODED RAW DATA IF NECESSARY */
	if (this_major_cnt > 1) {
	  if (fpout!=NULL) fwrite(obuff,sizeof(unsigned char),SAMPLES_PER_LINE,fpout);
	  else {printf("ERROR: Unable to write output to unopened file 1\n"); exit(1); }
          dump_all_headers(fp_hdr,this_major_cnt-1,last_major_sync_loc,&s);
	  for (i=0; i<SAMPLES_PER_LINE; i++) obuff[i] = 0;
	  optr = 0;
	}
    }
      
    /* frame #127 seems to be a sentinel...  I'm assuming this is NOT good data.
       So, make this the end of the dataset... */
    if (a.frame_no==127 && next_frame==127) {
         lock = 1;
	 error = 1;
     }

    if (DUMP_ALL_SYNCS) {
      double ftmp = this_sync; if (aligned==0) ftmp+=0.5;
      fprintf(fptmp,"%10.10i\t%10.1f\t%8.6i\t%3.3i\n",found_cnt,ftmp,major_cnt,a.frame_no);
    }

    if (headers < 10 && lock==1 && error==0) { 
        /*  This grabs every frame after the first frame_no 0 is hit,
	    i.e. it ignores the rest of the frame numbers, assuming that they 
	    are in sequence. This is done to avoid potential bit errors in 
	    the frame_no... */
	
	decode_headers(&r,buf,&headers);

	if (headers==10) {
	  decode_raw(&r,&s);
	
	  /* dump every single decoded time to this file */
          if (DUMP_TIMES==1) { fprintf(fp_time,"%li %lu\n",major_cnt,s.msec); }

	  /* dump every single thing that was decoded */
	  if (DUMP_ALL_HEADERS==1) { dump_all_headers(fp_all_hdr,major_cnt,major_sync_loc,&s); }

	  /* At the start of a data segment, dump a header */
	  // if (this_major_cnt == 1) {
	  //     int which = 0;
	  //     print_decoded_header(outheadername,major_cnt, major_sync_loc, &s, major_sync_num, which);
	  // }

	  /* every 1000 times we have all 9 aux headers filled, decode and display */
	  if ((major_cnt-1)%1000==0)  { 
	    display_decoded_header(major_cnt, major_sync_loc, &s, major_sync_num);
	    if (major_cnt==1) {
	      start_date.year = 1970+s.lsd_year;
	      start_date.jd = s.day_of_year;
	      dtmp = s.msec / 1000.0;
	      if (start_date.year != 1978) { 
	        major_cnt = 0; 
		this_major_cnt = 0;
		lock = 0;
		dataset_number--;
	        fclose(fpout);
	        fclose(fp_hdr);
	        fpout=NULL;
	        fp_hdr=NULL;
	        remove(outname);
	        remove(outheadername);
		optr = 0;
	        printf("Bad timing found at start of datatake... trying again...\n");
	      }
	    }
	  }
	}
    }  /* headers < 10 */
      
    /* because of bit errors throwing off the fill flag, we decode ALL frames as long as
       there is a lock...  It used to be this: if (a.fill_flag==0 && lock==1) {   */
      
    if (lock==1) {
	good_cnt++;
	int ret = decode_payload(buf,obuff,&optr);
	if (ret != 0)
	   printf("ERROR: Tried to write past end of obuff - discarding; line=%i a.frame_no=%i L3=%3.3i L2=%3.3i Last=%3.3i\n",
    	      major_cnt,a.frame_no,l3_frame,l2_frame,last_frame);
    }
    else if (a.fill_flag==1) {
	fill_cnt++;
	if (last_sync == last_bad) contiguous_miss++;
	else { 
	  max_consecutive_fills=(max_consecutive_fills>contiguous_miss) ? max_consecutive_fills : contiguous_miss;
	  contiguous_miss=0;
	}
	last_bad = this_sync;
        if (contiguous_miss > MAX_CONTIGUOUS_MISSES) end_of_dataset = 1;
    }
      
    if (lock==0) pre_cnt++;
    last_sync = this_sync;
    if (this_sync+(INT_FRAME_LEN*2) > file_length) {
	printf("MANUALLY predicted end of the file (%li bytes remain)!\n",file_length-this_sync);
	done=1;
    }
    
    if (error==1 && lock ==1) { /* we lost sync, dump data and try to establish it again */
      if (fpout!=NULL) {
         if (this_major_cnt > 10000) {
           fwrite(obuff,sizeof(unsigned char),SAMPLES_PER_LINE,fpout);
	   dump_all_headers(fp_hdr,this_major_cnt,major_sync_loc,&s);
	   fclose(fpout); 
	   fclose(fp_hdr);
	   fpout=NULL; 
	   fp_hdr=NULL;
	 
           /* decode_raw(&r,&s); */
	   display_decoded_header(major_cnt, major_sync_loc, &s, major_sync_num);
           // int which = this_major_cnt;
           // print_decoded_header(outheadername,major_cnt, major_sync_loc, &s, major_sync_num, which);

           printf("==========================================================================\n");
           printf("Lost Sync - Closed output file %s - dumped %i range lines\n",outname, this_major_cnt);
           printf("==========================================================================\n");
	 } else {  /* not enough data, throw it out */
	   dataset_number--;
	   fclose(fpout);
	   fclose(fp_hdr);
	   fpout=NULL;
	   fp_hdr=NULL;
	   remove(outname);
	   remove(outheadername);

           printf("==========================================================================\n");
           printf("Lost Sync - Destroyed output file %s - not enough lines %i\n",outname, this_major_cnt);
           printf("==========================================================================\n");
	 }
      } 
      
      for (i=0; i<SAMPLES_PER_LINE; i++) obuff[i] = 0;
      optr = 0;
      lock = 0;
      this_major_cnt=0;
      error = 0;
    }
   
    if (end_of_dataset == 1 && lock == 1) {
      printf("==========================================================================\n");
      printf("Consecutive fill - Closing output file %s - dumped %i range lines\n",outname,this_major_cnt-1);
      printf("==========================================================================\n");
      
      /*  Do not dump this last buffer, as it contains all fill data!
      if (fpout!=NULL) fwrite(obuff,sizeof(unsigned char),SAMPLES_PER_LINE,fpout);
      else {printf("ERROR: Unable to write output to unopened file 2\n"); exit(1); }
      */
      if (fpout!=NULL) { fclose(fpout); fclose(fp_hdr); fpout=NULL; fp_hdr=NULL; }
      for (i=0; i<SAMPLES_PER_LINE; i++) obuff[i] = 0;
      optr = 0;
      end_of_dataset=0;
      lock = 0;
      this_major_cnt=0;
    }
    
    /*******************************************************************************************
    This will break files into smaller pieces...
    Commented out because it was decided to decode entire swaths into single files
   
    if (this_major_cnt >= MAX_FILE_LENGTH && lock == 1) {
      fwrite(obuff,sizeof(unsigned char),SAMPLES_PER_LINE,fpout);
      fclose(fpout); 
      fpout=NULL; 
	 
      display_decoded_header(major_cnt, major_sync_loc, &s, major_sync_num);
      int which = this_major_cnt;
      print_decoded_header(outheadername,major_cnt, major_sync_loc, &s, major_sync_num, which);

      printf("==========================================================================\n");
      printf("FULL FILE - Closed output file %s - dumped %i range lines\n",outname, this_major_cnt);
      printf("==========================================================================\n");
      
      for (i=0; i<SAMPLES_PER_LINE; i++) obuff[i] = 0;
      optr = 0;
      lock = 0;
      this_major_cnt=0;
      error = 0;
    }
    ********************************************************************************************/
    
    
  }  /* while (!feof(fpin)&&!done&&!error) */

  /* Write out final line of data */
  if (optr != 0) {
    if (fpout!=NULL) {
      if (this_major_cnt > 10000) {
        fwrite(obuff,sizeof(unsigned char),SAMPLES_PER_LINE,fpout);
	dump_all_headers(fp_hdr,this_major_cnt,major_sync_loc,&s);
	fclose(fpout); 
	fclose(fp_hdr);
	fpout=NULL; 
	fp_hdr=NULL;
        decode_raw(&r,&s);
        // int which = this_major_cnt;
        // print_decoded_header(outheadername,major_cnt,major_sync_loc,&s,major_sync_num,which);
        printf("==========================================================================\n");
        printf("End of Data - Closed output file %s - dumped %i range lines\n",outname, this_major_cnt);
        printf("==========================================================================\n");
      } else {  /* not enough data, throw it out */
	fclose(fpout);
	fclose(fp_hdr);
	fpout=NULL;
	fp_hdr=NULL;
	remove(outname);
	remove(outheadername);
        printf("==========================================================================\n");
        printf("Lost Sync - Destroyed output file %s - not enough lines %i\n",outname, this_major_cnt);
        printf("==========================================================================\n");
      }
    }    
  }
  
  /* Display the last valid set of headers */
  decode_raw(&r,&s);
  display_decoded_header(major_cnt, major_sync_loc, &s, major_sync_num);

  printf("Found %i syncs in search; %i pre; %i good; %i fill (%i assumed non-fill)\n",found_cnt,pre_cnt,good_cnt,
  	fill_cnt,found_cnt-good_cnt-pre_cnt);
  printf("Found %i consecutive fill frames that were assumed to be good data\n",max_consecutive_fills); 
  printf("Fixed %i bad frame numbers (%i unfixed)\n",fixed_frame_nos, non_fixed_frame_nos);
  printf("Found %i partial lines\n",partial_line_cnt);
  printf("Wrote %i lines of output\n",major_cnt);

  fclose(fpin);

  if (DUMP_FIXED_FRAMES==1)     fclose(frame_file1);
  if (DUMP_NON_FIXED_FRAMES==1) fclose(frame_file2);
  if (DUMP_TIMES==1)            fclose(fp_time);
  if (DUMP_ALL_HEADERS==1)      fclose(fp_all_hdr);
  if (DUMP_ALL_SYNCS==1)        fclose(fptmp);
  
  exit(1);
}


void shift_buffer(unsigned char buf[], unsigned char first_val)
{
	unsigned char tbuf[INT_FRAME_LEN];
	int i;
	
	tbuf[0] = first_val << 4 | buf[0] >> 4;
	for (i=0; i<INT_FRAME_LEN-4; i++)  {   /* TAL - CHECK THIS -4 here!!! */
	  tbuf[i+1] = buf[i] << 4 | buf[i+1] >> 4;
	}
	for (i=0; i<INT_FRAME_LEN; i++) buf[i] = tbuf[i];
}
