

given: start_line, end_line

get_values(fphdr,hdr)  <-- first line

dwp_val[0]  = hdr->delay;
dwp_line[0] = 0;
dwp_cnt = 1;

/* go through the file and save all of the DWP changes */
for (i=start_line+1; i<end_line; i++) {
  get_values(fphdr,hdr)
  if (dwp_val[dwp_cnt-1] != hdr->delay) {
    dwp_val[dwp_cnt] = hdr->delay;
    dwp_line[dwp_cnt] = i;
    dwp_cnt++;
  }
}


if (dwp_cnt > 1) {  /* we have at least one DWP change */
  dwpfile = fopen(..., "w");

  dwp_min = 65;
  for (i=0;i<dwp_cnt; i++) { if (dwp_val[i] < dwp_min) dwp_min = dwp_val[i]; }

  if (dwp_min == dwp_val[0]) increasing = 1;  /* if increasing, first is zero - skip it */
  else increasing = 0;

  for (i=increasing; i<dwp_cnt; i++) {
      val = (dwp_val[i] - dwp_min)*216;
      fprintf(dwpfile,"%i %i\n",dwp_line[i],val);
  }
  
  fclose(dwpfile);
}
 
  
