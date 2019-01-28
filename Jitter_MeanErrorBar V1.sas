*options mlogic ;
/**********************************************************************************************
***********************************************************************************************;
Macro: %Jitter_MeanError
Created Date/Author: Jan 23, 2017/Chao Zhang
Current Version: V1
Working Environment: SAS 9.4 English version

Contact: Dr. Yuan Liu yliu31@emory.edu

Purpose:  To create a Jitter Mean with Standard Error for Longitudinal Data. This macro can be used
to create plots for batch of Y-vars(such as gene expression ) with same X-var.

Parameters: 
   DATA	        The name of the data set to be analyzed.
   X_VAR        Independent variable name (usually time or visit #); this is the variable to be jittered.                                                            
   Y_VAR	    Dependent variable name.                             
   GROUP_VAR    Grouping variable, id, or treatment allocation; must be numeric (optional).  
   XAXIS_VALUE	Usage values (example: (0 12 24 36) or (0 to 18 by 3)...; specify the ticket value of X axis (optional).  
   FINENAME	    File name for output table.  
   OUTPATH   	File path for output table to be stored.  
   DEBUG    	Set to T if running in debug mode (optional).  Work datasets will not be deleted in debug mode.  
                This is useful if you are editing the code or want to further manipulate the resulting data sets.  The default value is F.

******************************************************************************************************
******************************************************************************************************/
 
%macro Jitter_MeanError (data= , x_var=, y_var=, group_var=, xaxis_value=, outpath=, filename=,debug=F );

   %let debug = %UPCASE(&debug);
   ods rtf file="&OUTPATH.&FILENAME &SYSDATE..DOC.doc" startpage=NO;
   
   %let n=1;   
   %DO %UNTIL (%scan(&y_var, &n)= );
   %let y=%scan(&y_var, &n);

   
   ODS SELECT none ; 

   **To get label for y**;
   proc contents data=&data position out=vars; run;

   data vars;
         set vars;
         if label =" " then label = name;
		 if name="&y" then do;
		   label_y = label ;
         end;

		 if name="&x_var" then do;
		   label_x = label ;
         end;
	run;

	 
   proc sort data=vars; by descending label_y; run;

   data null;
     set vars;
	 call symputx("label_y", label_y);
	 if _n_=1;
    run;
  ** To get label for x**;
		 
   proc sort data=vars; by descending label_x; run;

   data null;
      set vars;
	  call symputx("label_x", label_x);
	  if _n_=1;
   run;
 
 
   %if &group_var ~= %STR() %then %do;

	  data wk;
	    set &data;
	    keep &y &x_var &group_var; 
	  run;
 
   proc sort data= wk;                                                                                                                    
      by &group_var  &x_var;                                                                                                                              
   run;  

   proc means data=wk noprint;                                                                                                           
      by  &group_var  &x_var;                                                                                                                               
      var &y;                                                                                                                            
      output out=meansout mean=mean stderr=stderr;                                                                                         
   run;     

  /* Calculate the upper and lower error bar values. */                                                                                       
   data reshape(drop=stderr);                                                                                                              
      set meansout;                                                                                                                        
      lower=mean - stderr;                                                                                                                 
      upper=mean + stderr;                                                                                                                 
   run;                                                                                                                                    
                                                                                                                                        
   proc freq data=reshape noprint;
       table &group_var  /out=freq_group; 
	   table &x_var  /out=freq_x;
   run;

   proc sql noprint;
      select &group_var into: count_G separated by ' ' from freq_group; 
   quit;
   /* Count number of freq in group */
   %let count_group = %sysfunc(countw(&count_G,' '));
   %put count_group :&count_group ;


   proc sql noprint;
      select &x_var into: count_X separated by ' ' from freq_X; 
   quit;
   /* Count number of freq in group */
   %let count_XVAR = %sysfunc(countw(&count_X,' '));
   %put count_XVAR :&count_XVAR ;

   data reshape; set reshape;
      adjust=5;
      if &count_XVAR > 9 then adjust=3;
      if &count_XVAR < 4 then adjust=10;
      i=1;
      if mod(&group_var,2)=0 then i=-1;
      &x_var = &x_var + (&count_group* i*(1/(adjust*&count_XVAR)));
    run;

    ods select all;
    title "Mean with Standard Error"; 
    proc sgplot data=reshape ;                                                                                                  
         scatter x=&x_var y=mean /group= &group_var yerrorlower=lower yerrorupper=upper                                                                 
                              markerattrs=( symbol=squareFilled size=12)  ;                                                                
         series x=&x_var y=mean / group = &group_var lineattrs=( pattern=solid) ;  

         xaxis label="&label_x" %if &xaxis_value ~= %STR() %then %do;
         values= &xaxis_value; %end;;
         yaxis label="&label_y" ; 
                                                       
    run; 
 %end;
    %else %do;
   
    data wk;
	    set &data;
	    keep &y &x_var ; 
	run;
 

    proc sort data= wk;                                                                                                                    
         by   &x_var;                                                                                                                              
    run;  

                                                                                                                                          
    proc means data=wk noprint;                                                                                                           
       by  &x_var;                                                                                                                               
       var &y;                                                                                                                            
       output out=meansout mean=mean stderr=stderr;                                                                                         
    run;     

    /* Calculate the upper and lower error bar values. */                                                                                       
    data reshape(drop=stderr);                                                                                                              
       set meansout;                                                                                                                        
       lower=mean - stderr;                                                                                                                 
       upper=mean + stderr;                                                                                                                 
    run; 
    ods select all;
    Title "Mean with Standard Error";
 
    proc sgplot data=reshape nonautolegend;                                                                                                  
         scatter x=&x_var y=mean / yerrorlower=lower yerrorupper=upper                                                                 
                            markerattrs=( symbol=squareFilled size=12)  ;                                                                
         series x=&x_var y=mean /  lineattrs=( pattern=solid) ;  
         xaxis label="&label_x" %if &xaxis_value ~= %STR() %then %do;
         values= &xaxis_value; %end;;
         yaxis label="&label_y" ; 
    run; 
%end;
   TITLE;


   %let n=%eval(&n+1);
   %end;
   ods rtf close;
/* If not in debug mode */

   %if &debug = F %then %do;
    PROC DATASETS lib=work noprint;
       delete vars  null Freq_group Freq_x meansout reshape ; 
    QUIT;
   %end;

%mend Jitter_MeanError;
