#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seasat.h"

/*************************************************************************************
Decode payload data into the output buffer
*************************************************************************************/
int decode_payload(unsigned char *buf, unsigned char *obuff, int *optr)
{
    int bufptr, j, which;

    bufptr = 1;  /*start at beginning of valid data */
    
    for (j=0;j<SAMPLES_PER_FRAME;j++)  {
      which = j%8;  /* which sample we are decoding */
      if (bufptr >= INT_FRAME_LEN+1) printf("ERROR: bufptr = %i past end of buffer!!! \n",bufptr);
      if ((*optr) >= SAMPLES_PER_LINE+SAMPLES_PER_FRAME) { return(1);   /* this is an error */ }
      
      switch (which) {
    	case 0:  obuff[*optr] = buf[bufptr] >> 3;
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 break;
    	case 1:  obuff[*optr] = ((buf[bufptr]&07)<<2) | (buf[bufptr+1]>>6);
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 bufptr++;
    		 break;
    	case 2:  obuff[*optr] = (buf[bufptr]&076)>>1;
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 break;
    	case 3:  obuff[*optr] = ((buf[bufptr]&01)<<4) | ((buf[bufptr+1]&0360)>>4);
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 bufptr++;
    		 break;
    	case 4:  obuff[*optr] = ((buf[bufptr]&017)<<1) | ((buf[bufptr+1]&0200)>>7);
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 bufptr++;
    		 break;
    	case 5:  obuff[*optr] = (buf[bufptr]&0174)>>2;
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 break;
    	case 6:  obuff[*optr] = ((buf[bufptr]&03)<<3) | ((buf[bufptr+1]&0340)>>5);
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 bufptr++;
    		 break;
    	case 7:  obuff[*optr] = buf[bufptr]&037;
    		 if (obuff[*optr]>31) {
    		    printf("ERROR: bad decode value %o (%i); case %i; obuff = %o\n",obuff[*optr],obuff[*optr],which,buf[bufptr]);
    		 }
    		 *optr+=1;
    		 bufptr++;
    		 break;
      }
    }
    return(0);
}

/*************************************************************************************
Convert raw buffer into raw aux headers
*************************************************************************************/
void decode_headers(SEASAT_raw_header *r, unsigned char *buf, int *headers)
{
	switch (*headers) {
	  case 0:  
	    	r->aux0->lsd_year = buf[0] >> 4;
		r->aux0->station  = buf[0] & 017;
		*headers+=1; 
		/* display_aux(&r,0); */
		break;   
	  case 1:  
	  	r->aux1->msec0 = buf[0];
		*headers+=1; 
		/* display_aux(&r,1); */
	  	break;      
	  case 2:  
	        r->aux2->msec8 = buf[0];
		*headers+=1;
		/* display_aux(&r,2); */
	  	break;      
	  case 3:
	        r->aux3->msec16 = buf[0];
		*headers+=1; 
		/* display_aux(&r,3); */
	  	break;      
	  case 4:  
	        r->aux4->day_of_year0 = buf[0] >> 3;
		r->aux4->msec24 = buf[0] & 007;
		*headers+=1; 
		/* display_aux(&r,4); */
	  	break;      
	  case 5:  
	        r->aux5->clock_drift0 = buf[0] >> 4;
		r->aux5->day_of_year5 = buf[0] & 017;
		*headers+=1; 
		/* display_aux(&r,5); */
	  	break;      
	  case 6:
	        r->aux6->clock_drift4 = buf[0];
		*headers+=1; 
		/* display_aux(&r,6); */
	  	break;      
	  case 7:
	        r->aux7->scan_indicator = buf[0] >> 7;
		r->aux7->bits_per_sample = (buf[0] & 0160) >> 4;
		r->aux7->mfr_lock = (buf[0] & 010) >> 3;
		r->aux7->prf_code = buf[0] & 007;
		*headers+=1; 
		/* display_aux(&r,7); */
	  	break;      
	  case 8:  
	        r->aux8->delay_10 = buf[0] >> 4;
		r->aux8->delay_1  = buf[0] & 017;
		*headers+=1; 
		/* display_aux(&r,8); */
	  	break;      
	  case 9:  
	        r->aux9->scu_bit = (buf[0] & 200) >> 7;
		r->aux9->sdf_bit = (buf[0] & 100) >> 6;
		r->aux9->adc_bit = (buf[0] & 040) >> 5;
		r->aux9->time_gate_bit = (buf[0] & 020) >> 4;
		r->aux9->local_prf_bit = (buf[0] & 010) >> 3;
		r->aux9->auto_prf_bit = (buf[0] & 004) >> 2;
		r->aux9->prf_lock_bit = (buf[0] & 002) >> 1;
		r->aux9->local_delay_bit = buf[0] & 001;
		*headers += 1;
		/* display_aux(&r,9);  */
	  	break;      
	}
}



/*************************************************************************************
  Populate the decoded structure from raw aux values
*************************************************************************************/
void decode_raw(SEASAT_raw_header *r, SEASAT_header *s) {
  s->lsd_year = r->aux0->lsd_year;
  s->station_code = r->aux0->station;
  s->msec= r->aux1->msec0 |
  	   r->aux2->msec8 << 8 |
	   r->aux3->msec16 << 16 | 
	   r->aux4->msec24 << 24;
  s->day_of_year = r->aux4->day_of_year0 |
  		   r->aux5->day_of_year5 << 5;
  s->clock_drift = r->aux5->clock_drift0 |
  		   r->aux6->clock_drift4 << 4;
  s->no_scan_indicator_bit = r->aux7->scan_indicator;
  s->bits_per_sample = r->aux7->bits_per_sample;
  s->mfr_lock_bit = r->aux7->mfr_lock;
  s->prf_rate_code = r->aux7->prf_code;
  s->delay = 10*r->aux8->delay_10 + r->aux8->delay_1;
  s->scu_bit = r->aux9->scu_bit;
  s->sdf_bit = r->aux9->sdf_bit;
  s->adc_bit = r->aux9->adc_bit;
  s->time_gate_bit = r->aux9->time_gate_bit;
  s->local_prf_bit = r->aux9->local_prf_bit;
  s->auto_prf_bit  = r->aux9->auto_prf_bit;
  s->prf_lock_bit  = r->aux9->prf_lock_bit;
  s->local_delay_bit = r->aux9->local_delay_bit;
}

/*************************************************************************************
  Display the decoded structure
*************************************************************************************/
void display_decoded_header(int major_cnt, long int this_sync_loc, SEASAT_header *s, int major_sync_num) {
  printf("\n------------------ HEADER %i INFORMATION (SYNC # %i)------------------\n",major_cnt,major_sync_num);
  printf("Sync offset for 0th minor frame %li\n",this_sync_loc);
  printf("Least significant digit of year  %i\n",s->lsd_year);
  printf("Station is");
  if (s->station_code==5) printf(" Alaska ");
  else if (s->station_code==6) printf(" Goldstone ");
  else if (s->station_code==7) printf(" Merrit Island ");
  else if (s->station_code==9) printf(" Oak Hangar ");
  else if (s->station_code==10) printf(" Shoe Code ");
  else printf(" UNKNOWN ");
  printf("(code is %i)\n",s->station_code);
  printf("Time of day is %lu\n",s->msec);  /* Go back and make this a time TOM! */
  printf("Day of the year %i\n",s->day_of_year); /* make into date */
  printf("Sample Clock drift %i\n",s->clock_drift);
  printf("No scan indicator bit  %i\n",s->no_scan_indicator_bit);
  printf("Bits per sample %i\n",s->bits_per_sample);
  printf("MFR lock bit %i\n",s->mfr_lock_bit);
  printf("PRF rate");
  if (s->prf_rate_code==1) printf(" 1464 ");
  else if (s->prf_rate_code==2) printf(" 1540 ");
  else if (s->prf_rate_code==3) printf(" 1581 ");
  else if (s->prf_rate_code==4) printf(" 1647 ");
  else printf(" INVALID CODE ");
  printf("(code is %i)\n",s->prf_rate_code);
  printf("Delay from leading edge of PRF to start of digitization %i\n",s->delay);
  printf("SCU bit %i\n",s->scu_bit);
  printf("SDF bit %i\n",s->sdf_bit);
  printf("ADC bit %i\n",s->adc_bit);
  printf("Time Gate bit %i\n",s->time_gate_bit);
  printf("Local PRF bit %i\n",s->local_prf_bit);
  printf("Auto PRF bit %i\n",s->auto_prf_bit);
  printf("PRF Lock bit %i\n",s->prf_lock_bit);
  printf("Local Delay bit %i\n",s->local_delay_bit);
  printf("------------------ END HEADER INFORMATION ------------------\n\n");

}

/*************************************************************************************
  Print the decoded structure to file
*************************************************************************************/
void print_decoded_header(char *outheadername,int major_cnt,long int this_sync_loc,
			SEASAT_header *s,int major_sync_num, int which) 
{
  FILE *fp;
  
  /* switch==0 means start a new file */
  if (which==0) fp = fopen(outheadername,"w");
  else fp = fopen(outheadername,"a");
  if (fp == NULL) {printf("Error opening header file %s\n",outheadername); exit(1);}

  fprintf(fp,"\n------------------ HEADER: %i INFORMATION (SYNC # %i)------------------\n",major_cnt,major_sync_num);
  fprintf(fp,"Sync offset for 0th minor frame: %li\n",this_sync_loc);
  fprintf(fp,"Least significant digit of year: %i\n",s->lsd_year);
  fprintf(fp,"Station is:");
  if (s->station_code==5) fprintf(fp," Alaska ");
  else if (s->station_code==6) fprintf(fp," Goldstone ");
  else if (s->station_code==7) fprintf(fp," Merrit Island ");
  else if (s->station_code==9) fprintf(fp," Oak Hangar ");
  else if (s->station_code==10) fprintf(fp," Shoe Code ");
  else fprintf(fp," UNKNOWN ");
  fprintf(fp,"(code is %i)\n",s->station_code);
  fprintf(fp,"Time of day is: %lu\n",s->msec);  /* Go back and make this a time TOM! */
  fprintf(fp,"Day of the year: %i\n",s->day_of_year); /* make into date */
  fprintf(fp,"Sample Clock drift: %i\n",s->clock_drift);
  fprintf(fp,"No scan indicator bit: %i\n",s->no_scan_indicator_bit);
  fprintf(fp,"Bits per sample: %i\n",s->bits_per_sample);
  fprintf(fp,"MFR lock bit: %i\n",s->mfr_lock_bit);
  fprintf(fp,"PRF rate:");
  if (s->prf_rate_code==1) fprintf(fp," 1464 ");
  else if (s->prf_rate_code==2) fprintf(fp," 1540 ");
  else if (s->prf_rate_code==3) fprintf(fp," 1581 ");
  else if (s->prf_rate_code==4) fprintf(fp," 1647 ");
  else fprintf(fp," INVALID CODE ");
  fprintf(fp,"(code is %i)\n",s->prf_rate_code);
  fprintf(fp,"Delay from leading edge of PRF to start of digitization: %i\n",s->delay);
  fprintf(fp,"SCU bit: %i\n",s->scu_bit);
  fprintf(fp,"SDF bit: %i\n",s->sdf_bit);
  fprintf(fp,"ADC bit: %i\n",s->adc_bit);
  fprintf(fp,"Time Gate bit: %i\n",s->time_gate_bit);
  fprintf(fp,"Local PRF bit: %i\n",s->local_prf_bit);
  fprintf(fp,"Auto PRF bit: %i\n",s->auto_prf_bit);
  fprintf(fp,"PRF Lock bit: %i\n",s->prf_lock_bit);
  fprintf(fp,"Local Delay bit: %i\n",s->local_delay_bit);
  fprintf(fp,"------------------ END HEADER INFORMATION ------------------\n\n");
  
  if (which!=0) fprintf(fp,"Number of lines: %i\n",which);
  fclose(fp);
}


/*************************************************************************************
  Display raw aux data
*************************************************************************************/
void display_aux(SEASAT_raw_header *r, int field){
 switch (field) {
  case 0: printf("AUX 0: last digit of year %i, station code %i\n",
   	  r->aux0->lsd_year,r->aux0->station); break;
  case 1: printf("AUX 1: msec (lower 8 bits) %i\n",r->aux1->msec0); break;
  case 2: printf("AUX 2: msec (starting at bit 8) %i\n",r->aux2->msec8); break;
  case 3: printf("AUX 3: msec (starting at bit 16) %i\n",r->aux3->msec16); 
   	  break;
  case 4: printf("AUX 4: day of year (lower 5 bits) %i, msec (starting at bit 24) %i\n",
   	  r->aux4->day_of_year0,r->aux4->msec24); break;
  case 5: printf("AUX 5: clock drift (lower 4 bits) %i, day of year (upper 4 bits) %i\n",
  	  r->aux5->clock_drift0,r->aux5->day_of_year5); break; 
  case 6: printf("AUX 6: clock drift (upper 8 bits) %i\n",r->aux6->clock_drift4); break;
  case 7: printf("AUX 7: scan_id %i, bits per sample %i, mfr_lock %i, prf_code %i\n",
  	  r->aux7->scan_indicator, r->aux7->bits_per_sample,
	  r->aux7->mfr_lock, r->aux7->prf_code); break;
  case 8: printf("AUX 8: delay10 %i, delay1 %i\n",r->aux8->delay_10, 
  	  r->aux8->delay_1); break;
  case 9: printf("AUX 9: Misc bits= scu:%i sdf:%i adc:%i time_gate:%i\n",
 	  r->aux9->scu_bit,r->aux9->sdf_bit,r->aux9->adc_bit,r->aux9->time_gate_bit);
  	  printf("AUX 9: Misc bits= local:%i auto:%i lock:%i delay:%i\n",
	  r->aux9->local_prf_bit,r->aux9->auto_prf_bit,r->aux9->prf_lock_bit,
	  r->aux9->local_delay_bit); break;
 }
}

 
/*************************************************************************************
  Dump all decoded header values to files
*************************************************************************************/
void dump_all_headers(FILE *fp_all_hdr,int major_cnt,long int major_sync_loc,SEASAT_header *s)
{
  fprintf(fp_all_hdr,"%i %li %i %i %i %li %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n",
    major_cnt,
    major_sync_loc,
    s->station_code,
    s->lsd_year,
    s->day_of_year,
    s->msec,
    s->clock_drift,
    s->no_scan_indicator_bit,
    s->bits_per_sample,
    s->mfr_lock_bit,
    s->prf_rate_code,
    s->delay,
    s->scu_bit,
    s->sdf_bit,
    s->adc_bit,
    s->time_gate_bit,
    s->local_prf_bit,
    s->auto_prf_bit,
    s->prf_lock_bit,
    s->local_delay_bit);
}

