#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seasat.h"

unsigned char SYNC1='\371';
unsigned char SYNC2='\250';
unsigned char SYNC3='\355';

/*-------------------------------------------------------------------------
 Searches for a sync code starting at current location of the file pointer.
 Returns 1 if found; 0 if not.
 
 This routine may have problems at the end of the file...  I never check
 for end of file conditions...
 -------------------------------------------------------------------------*/
int find_sync(FILE *fpin)
{
	unsigned char tmp1, tmp2, tmp3;
	int found;
	
	found = 0;
	while (found == 0) {
	  while(tmp1 != SYNC1) if (fread(&tmp1,sizeof(unsigned char),1,fpin)!=1) return(-1);
	  if (fread(&tmp2,sizeof(unsigned char),1,fpin)!=1) return(-1);
	  if (tmp2 != SYNC2) continue;
	  if (fread(&tmp3,sizeof(unsigned char),1,fpin)!=1) return(-1);
	  if (tmp3 != SYNC3) continue;
	  printf("Found sync code at %li\n",ftell(fpin)-3);
	  found=1;
	}
	return(found);
}

/*-------------------------------------------------------------------------
   Looks for a sync code at the current location
   Returns 1 if found, 0 if not found, -1 if EOF
   Allows for BIT_ERRORS # of errors when looking for a match.
 -------------------------------------------------------------------------*/
int find_one_sync(FILE *fpin)
{
	unsigned char tmp1, tmp2, tmp3;
	int found;
	int be;
	
	found = 0;
	be = 0;
	if (fread(&tmp1,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if( fread(&tmp2,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if (fread(&tmp3,sizeof(unsigned char),1,fpin)!=1) return(-1);
  	if (tmp1 == SYNC1 && tmp2 == SYNC2 && tmp3 == SYNC3) found = 1;
	else {
	  if (tmp1 != SYNC1) be = bit_errors(SYNC1,tmp1);
	  if (tmp2 != SYNC2) be += bit_errors(SYNC2,tmp2);
	  if (tmp3 != SYNC3) be += bit_errors(SYNC3,tmp3);
	  if (be<=BIT_ERRORS) found = 1;
          /* else printf("Error sync code at %.8i: %o %o %o (%i bit errors)\n",
	     ftell(fpin)-3,tmp1,tmp2,tmp3,be); */
	}
	return(found);
}
	
/* This routine finds the number of bits different between ref(erence) and pat(tern) */
int bit_errors(unsigned char ref, unsigned char pat)
{	
    int be, len, test, i;

    len =sizeof(ref);
    be = 0;
    test = ref ^ pat;
    for (i=0; i<len*8; i++)
      {
        int a;
	a = test & 1;
	be = be + a;
 	test = test >> 1;
      }
    return(be);
}

/*-------------------------------------------------------------------------
   Looks for a sync code at the current location
   Returns 1 if found, 0 if not found, -1 if EOF
   Allows for BIT_ERRORS # of errors when looking for a match.
 -------------------------------------------------------------------------*/
int find_sync_no_advance(FILE *fpin, int *bers)
{
	unsigned char tmp1, tmp2, tmp3;
	int found;
	int be;
	long offset;
	
	found = 0;
	be = 0;
	if (fread(&tmp1,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if( fread(&tmp2,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if (fread(&tmp3,sizeof(unsigned char),1,fpin)!=1) return(-1);
  	if (tmp1 == SYNC1 && tmp2 == SYNC2 && tmp3 == SYNC3) found = 1;
	else {
	  if (tmp1 != SYNC1) be = bit_errors(SYNC1,tmp1);
	  if (tmp2 != SYNC2) be += bit_errors(SYNC2,tmp2);
	  if (tmp3 != SYNC3) be += bit_errors(SYNC3,tmp3);
	  if (be<=BIT_ERRORS) found = 1;
          /* else printf("Error sync code at %.8i: %o %o %o (%i bit errors)\n",
	     ftell(fpin)-3,tmp1,tmp2,tmp3,be); */
	}
	
        offset = ftell(fpin);
	offset = offset - 3;
        fseek(fpin,offset,SEEK_SET);
	*bers = be;
	return(found);
}

/*-------------------------------------------------------------------------
   Looks for a sync code at the current location
   Returns 1 if found, 0 if not found, -1 if EOF
   Allows for BIT_ERRORS # of errors when looking for a match.
 -------------------------------------------------------------------------*/
int find_unaligned_sync_no_advance(FILE *fpin, int *bers)
{
	unsigned char rtmp1, rtmp2, rtmp3, rtmp4;
	unsigned char tmp1, tmp2, tmp3;
	int found;
	int be;
	long offset;
	
	found = 0;
	be = 0;
	if (fread(&rtmp1,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if (fread(&rtmp2,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if (fread(&rtmp3,sizeof(unsigned char),1,fpin)!=1) return(-1);
	if (fread(&rtmp4,sizeof(unsigned char),1,fpin)!=1) return(-1);
	
	/* Shift the bits and combine into regular bytes for comparison */
	tmp1 = rtmp1 << 4 | rtmp2 >> 4;
	tmp2 = rtmp2 << 4 | rtmp3 >> 4;
	tmp3 = rtmp3 << 4 | rtmp4 >> 4;
		
  	if (tmp1 == SYNC1 && tmp2 == SYNC2 && tmp3 == SYNC3) found = 2;
	else {
	  if (tmp1 != SYNC1) be = bit_errors(SYNC1,tmp1);
	  if (tmp2 != SYNC2) be += bit_errors(SYNC2,tmp2);
	  if (tmp3 != SYNC3) be += bit_errors(SYNC3,tmp3);
	  if (be<=BIT_ERRORS) found = 2;
          /* else printf("Error unaligned sync code at %.8i: %o %o %o (%i bit errors)\n",
	     ftell(fpin)-4,tmp1,tmp2,tmp3,be); */
	}
	
        offset = ftell(fpin);
	offset = offset - 4;
        fseek(fpin,offset,SEEK_SET);
	*bers = be;
	return(found);
}

int get_next_frameno(FILE *fpin, int aligned)
  {
      long save_loc, this_loc;
      unsigned char tmpc, rtmp1, rtmp2;
      int  frame_no;
      
      save_loc = ftell(fpin);
      this_loc = save_loc + 3;
      fseek(fpin,this_loc,SEEK_SET);

      /* If we had (last) an unaligned sync code, do this to get the frame number */
      if (aligned==0) {
	if (fread(&tmpc,sizeof(unsigned char),1,fpin)!=1) {printf("ERROR reading frame no\n"); exit(1);}
	frame_no = tmpc & 0177;
      }

      /* If we had (last) an aligned sync code, do this to get the frame number */  	
      if (aligned==1 || aligned == 2) {
	if (fread(&rtmp1,sizeof(unsigned char),1,fpin)!=1) {printf("ERROR reading frame no\n"); exit(1);}
	if (fread(&rtmp2,sizeof(unsigned char),1,fpin)!=1) {printf("ERROR reading frame no\n"); exit(1);}
	tmpc = rtmp1 << 4 | rtmp2 >> 4;
	frame_no = tmpc & 0177;
      }
      
      fseek(fpin,save_loc,SEEK_SET);
      return(frame_no);
  }



























