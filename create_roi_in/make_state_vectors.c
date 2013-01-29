
  julian_date start_date;
  julian_date end_date;
  ymd_date    s_ymd;
  ymd_date    e_ymd;
  hms_time    start_time;
  hms_time    end_time;


 start_date.year = 1970+s.lsd_year;
 start_date.jd = s.day_of_year;
 dtmp = s.msec / 1000.0;
 date_sec2hms(dtmp,&start_time);

  end_date.year = 1970+s.lsd_year;
  end_date.jd = s.day_of_year;
  dtmp = s.msec / 1000.0;
  date_sec2hms(dtmp,&end_time);

  date_jd2ymd(&start_date,&s_ymd);
  date_jd2ymd(&end_date,&e_ymd);
  
  printf("Start of dataset (mm/dd/yyyy) %.2i/%.2i/%i %.2i:%.2i:%6.3lf\n",s_ymd.month,s_ymd.day,s_ymd.year,
  	start_time.hour,start_time.min,start_time.sec);

  printf("End   of dataset (mm/dd/yyyy) %.2i/%.2i/%i %.2i:%.2i:%6.3lf\n",e_ymd.month,e_ymd.day,e_ymd.year,
  	end_time.hour,end_time.min,end_time.sec);

  printf("\n\nCalculating state vectors\n");
  create_input_tle_file(start_date,start_time,"tle1.txt");
  propagate_state_vector("tle1.txt"); 
  printf("\n\nConverting state vectors from ECI to ECEF\n");
  fix_state_vectors(start_date.year,start_date.jd,start_time.hour,start_time.min,start_time.sec);
  
