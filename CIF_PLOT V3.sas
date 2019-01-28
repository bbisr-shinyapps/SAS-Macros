/**********************************************************************************************
***********************************************************************************************;
Macro Name: CIF_PLOT
Created Date/Author/Contact: Oct 2016/ Yaqi Jia, Chao Zhang
Last Update Date/Person: Oct, 2016/Yuan Liu
Other Signficant Contributor:
Current Version: V3
Working Environment: SAS 9.4

Contact Information: Dr. Yuan Liu/ yliu31@emory.edu

Purpose:  To create a cumulative incidence plot for competing risk analysis based on Fine&Gray model (1999). 
Generate time specific event rate with 95% confidence interval and Fine and Gray p value for comparison in term of 
event of interest. 

Note: All calculation is based on PROC PHREG in SAS 9.4 with EVENTCODE option. See SAS help manual for more 
technique details. 

Parameters: 

DSN            The name of the data set to be analyzed.   
 
GRPLIST        The variable list that defines the groups for comparison (optional) separated by space.

TIME_EVENT     Time to event outcome variable. 

CENSOR         Name of censoring indicator variable.  Values of 0 indicate censored observations.

EVENTCODE      Values in CENSOR parameter that correspond to the event of interest. The default value is 1.

YAXISVALUE     Specify the ticket value of y axis. 

XAXISVALUE	   Specify the ticket value of X axis.; 

TITLE1         First title for the plot  (optional).  

TITLE2         Sub title for the plot  (optional).
 
TIMELIST       List of time points separated by spaces to report survival estimates and 95% CI 
               (optional).  

UNITS          Units of the time variable, i.e. days, months, etc. The default value is none.
OUTPATH        Path for output table to be stored.
FNAME          File name for output table.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted in
               debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

******************************************************************************************************/
%macro CIF_PLOT(dsn=,grplist=,time_event=,censor=,eventcode=,title1=,title2=,xaxisvalue=,yaxisvalue=, timelist=,units=,outpath=,filename=,debug=F);

    
   /* Count number of time estimates requested */
   %let time_cnt = %sysfunc(countw(&timelist,' '));

/*count number of &grplist*/
    %let var_cnt = %sysfunc(countw(&grplist,' '));

ods rtf file="&OUTPATH.&FILENAME &SYSDATE..DOC.doc" startpage=NO;

%if &var_cnt = 0 %then %do;

		ODS SELECT NONE;
           proc phreg data=&dsn simple plots (overlay=stratum )=cif ;
                model &time_event*&censor(0) =  /eventcode=&eventcode;
                baseline out=plot_ cif=_all_;
           run;
		ODS SELECT ALL;

     %if &title1 ~= %STR() %then %do; TITLE1 &title1; %end;
     %if &title1 ~= %STR() %then %do; TITLE2 &title2; %end;

      ods select all;
      proc sgplot data= plot_;
	      %if &yaxisvalue ~= %str() %then %do; yaxis values=(&yaxisvalue);%end;
          %if &xaxisvalue ~= %str() %then %do; xaxis values= (&xaxisvalue);%end;
		step y=CIF x=&time_event ;
     run;
		
        %if &time_cnt > 0 %then %do;
		    DATA _est;
		        set plot_;
		        %do i = 1 %to &time_cnt;
		          %let timept = %SCAN(&timelist,&i,' ');
		         
		         /* Create indicator if time is less than or equal to time of interest */
		          if &time_event <= &timept then _LT&i = 1;
		         else _LT&i = 0;
		       %end;
		     RUN;

     PROC SQL;
        create table _est2
        as select *,  max(&time_event) as max_time, max(&time_event*_LT1) as _max1 
        %do i = 2 %to &time_cnt;  , max(&time_event*_LT&i) as _max&i %end;
        from _est;
        /*group by &censor;*/
     QUIT;

     DATA _est3;
        set _est2;
       /* Subset on time points of interest */
       %do i = 1 %to &time_cnt;
          %let timept = %SCAN(&timelist,&i,' ');

          if &time_event = _max&i then do;
            _time = &timept;
            output;
          end;
       %end;

     RUN;

     DATA _est3;
        set _est3;
       /* Convert to character */
       cif_c = PUT(cif,percent8.2);
       lcl_c = PUT(LowerCIF,percent8.2);
       ucl_c = PUT(upperCIF,percent8.2);

       /* Make sure estimates are not for times beyond the length of the data */
       if _time > max_time then do;
            cif_c = 'NA';
            Lcl_c = 'NA';
            Ucl_c = 'NA';
       end;

       /* Combine estimate and CI */
       rate = TRIM(cif_c) || " (" || TRIM(LEFT(LCL_c)) || ", " || TRIM(LEFT(UCL_c)) || ")";

     RUN;

     PROC SORT DATA = _est3; by  _time; RUN;

     proc report nowd data=_est3 style(report)={rules=GROUPS}  
            style(header)={BACKGROUNDCOLOR=none} LS=256 STYLE(COLUMN) = {JUST = C};
         col  _time rate;
         %if &units ~= %str() %then %do;
         LABEL _time = "Time (&units)"
		 %end;
		 %else %do; LABEL _time="Time"  %end;
            rate = 'CIF Estimate (95% CI)'; 
         RUN;

	%end; /*%if &time_cnt > 0 %then %do;*/
%end; /*End of %if &var_cnt = 0 %then %do;*/

%if &var_cnt > 0 %then %do;

    %do n = 1 %to &var_cnt;
		%let grp=%scan(&grplist,&n);
		ODS SELECT none ; 
        proc freq data=&dsn;
        table &grp/out=risk1;
        run;

		ODS SELECT NONE;
		ODS OUTPUT "Maximum Likelihood Estimates of Model Parameters" = _mlec  ClassLevelFreq =freq  type3=type3 ; 
           proc phreg data=&dsn simple plots (overlay=stratum )=cif ;
                class &grp / order=internal ref=first param=glm;
                model &time_event*&censor(0) = &grp  /eventcode=&eventcode;
                baseline covariates=risk1 out=plot_ cif=_all_/ rowid=&grp;
           run;
		 ODS SELECT ALL;

          Proc sql noprint;
             select ProbChiSq format=pvalue6.4 into:pval1
             from type3;
          quit;
     

        %if &title1 ~= %STR() %then %do; TITLE1 &title1; %end;
        %if &title1 ~= %STR() %then %do; TITLE2 &title2; %end;

        ods select all;
 
         proc sgplot data= plot_;
	      %if &yaxisvalue ~= %str() %then %do; yaxis values=(&yaxisvalue) label = "&ylab";%end;
          %if &xaxisvalue ~= %str() %then %do; xaxis values= (&xaxisvalue) label = "&xlab";%end;
			step y=CIF x=&time_event /group=&grp;
          	inset ("Fine-Gray p-value="="%TRIM(&pval1)")/position=topright border LABELALIGN=LEFT;
         run;

        %if &time_cnt > 0 %then %do;
		    DATA _est;
		        set plot_;
		        %do i = 1 %to &time_cnt;
		          %let timept = %SCAN(&timelist,&i,' ');
		         
		         /* Create indicator if time is less than or equal to time of interest */
		          if &time_event <= &timept then _LT&i = 1;
		         else _LT&i = 0;
		       %end;
		     RUN;

     PROC SQL;
        create table _est2
        as select *,  max(&time_event) as max_time, max(&time_event*_LT1) as _max1 
        %do i = 2 %to &time_cnt;  , max(&time_event*_LT&i) as _max&i %end;
        from _est;
        /*group by &censor;*/
     QUIT;

     DATA _est3;
        set _est2;
       /* Subset on time points of interest */
       %do i = 1 %to &time_cnt;
          %let timept = %SCAN(&timelist,&i,' ');

          if &time_event = _max&i then do;
            _time = &timept;
            output;
          end;
       %end;

     RUN;

     DATA _est3;
        set _est3;
       /* Convert to character */
       cif_c = PUT(cif,percent8.2);
       lcl_c = PUT(LowerCIF,percent8.2);
       ucl_c = PUT(upperCIF,percent8.2);

       /* Make sure estimates are not for times beyond the length of the data */
       if _time > max_time then do;
            cif_c = 'NA';
            Lcl_c = 'NA';
            Ucl_c = 'NA';
       end;

       /* Combine estimate and CI */
       rate = TRIM(cif_c) || " (" || TRIM(LEFT(LCL_c)) || ", " || TRIM(LEFT(UCL_c)) || ")";

     RUN;

     PROC SORT DATA = _est3; by &grp _time; RUN;

     proc report nowd data=_est3 style(report)={rules=GROUPS}  
            style(header)={BACKGROUNDCOLOR=none} LS=256 STYLE(COLUMN) = {JUST = C};
         col &grp _time rate;
         define &grp/order style={just=l};

         %if &units ~= %str() %then %do;
         LABEL _time = "Time (&units)"
		 %end;
		 %else %do; LABEL _time="Time"  %end;
            rate = 'CIF Estimate (95% CI)'; 
         RUN;


		%end; /*%if &time_cnt > 0 %then %do;*/
	%end;/* %do n = 1 %to &cvar_cnt;*/
%end; /*%if &var_cnt > 0 %then %do;*/
 
ods rtf close;
/* If not in debug mode */
   %if &debug = F %then %do;
      /* Delete intermediate datasets */
      PROC DATASETS lib=work noprint;
        delete Freq plot_ risk1 type3 
         %if &time_cnt > 0 %then %do; _est _est2 _est3 ;%end; 
      QUIT;
   %end;

%mend CIF_PLOT ;
