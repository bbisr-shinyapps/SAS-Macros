/**********************************************************************************************
Macro: LINREG_SEL
Created Date/Author: Oct. 3, 2013/Dana Nickleach
Last Update Date/Person: 
Current Version: V1
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To conduct backward selection on a general linear model using the maximum possible 
sample size at each stage of the selection process instead of restricting to the sample size 
from the first step as SAS does when using their selection methods.  Optionally, a table of the
resulting model can be generated.  

Notes: The model is run using PROC GLMSELECT and PROC GLM.  The final list of variables selected 
will be written to the log.  Additionally, two global macro variables, _finalvar and _finalcvar 
will be created containing the list of all variables and categorical variables selected, 
respectively.  If you are requesting a table with the model results then the macro 
“MUTLIPLE_LINREG V2” or later is also required.  Interactions can be included if you are not 
requesting a table (REPORT=F).  Interactions will be removed from the model with respect to 
model hierarchy.  However, the macro cannot produce a correct report table if the model contains 
interactions.  

Parameters: 

DSN            The name of the data set to be analyzed.
OUTCOME        The name of the continuous outcome variable.  
VAR            List of variables to include in the model separated by spaces.  Each variable 
               name or interaction term must not be more than 20 characters long.                                           
CVAR           List of categorical variables to include in the model separated by spaces.  
               These should also appear in the VAR parameter.
ORDER          The order in which to sort the levels of the categorical variables: DATA, 
               FORMATTED, FREQ, or INTERNAL.  The default value is INTERNAL.
INC            Number of variables to include in the model.  The first n variables in the 
               var parameter will be included in every model.  The default value is 0.  
SLSTAY         The significance level for removing variables from the model.  The default value
               is .05.   
WEIGHT         Variable to use in the weight statement.  The macro will not normalize weights to
               the original sample size.  Weights will need to be normalized in the data
               preparation stage.  Leave it blank if not using weights.
REPORT         Set this to T if you want a table of the resulting model generated.  The default 
               value is F.   
TYPE3          Set to F to suppress type III p-values from being reported in the table.  The 
               default value is T.  This only has an effect if REPORT = T.
DEBUG          Set to T if running in debug mode.  Work datasets will not be deleted in debug 
               mode.  This is useful if you are editing the code or want to further manipulate 
               the resulting data sets.  The default value is F.
OUTPATH        File path for output table to be stored. 
FILENAME       File name for output table.         

***********************************************************************************************
* For more details, please see the related documentation
**********************************************************************************************/

%macro linreg_sel(dsn=,outcome=,var=,cvar=,order=internal,inc=0,slstay=.05,weight=,report=F,
     type3=T,outpath=,filename=,debug=F);

   /* Macros for final variable lists */
   %global _finalvar _finalcvar;

   %local FILENAME FOOTTEXT WEIGHT DSN REPORT OUTPATH CVAR_CNT CVAR I REMOVE OUTCOME CONTINUE
     __MACRO_ERR EVENT VAR DESC INC TYPE3 REMLIST DEBUG SLSTAY VAR_CNT;

   /* Upper case T/F */
   %let report = %UPCASE(&report);
   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);
   %let var = %UPCASE(&var);
   %let cvar = %UPCASE(&cvar);

   /* Count number of variables */
   %let var_cnt = %sysfunc(countw(&var,' '));
   /* Count number of class variables */
   %let cvar_cnt = %sysfunc(countw(&cvar,' '));       

   /* Initialize error variable */
   %LET __Macro_Err = 0;

   /* Check for variable names that are too long */
   %DO i = 1 %to &var_cnt; 
      %IF %LENGTH(%SCAN(&var, &i, %STR( ))) > 20 %then %do;
         %put ERROR: Variable name %SCAN(&var, &i, %STR( )) must be less than 21 characters.; 
           %let __Macro_Err=1;
       %END;

      /* Check for interaction and report request */
       %if &report = T %then %do;
          %IF %INDEX(%SCAN(&var, &i, %STR( )),*) > 0 %then %do;
            %put ERROR: A report cannot be created if interaction terms are in the model.  Set REPORT=F.; 
              %let __Macro_Err=1;
          %END;
       %end;
   %END;

   /* If there is an error in the parameters supplied then exit */
   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Intialize */
   %let continue = 1;
   /* List of variables removed */
   %let remlist =;

   /* Run one step of backwards selection at a time and stop once there are no more variables*/
   /*  that should be removed. */
   %do %while(&continue=1);

   
      ODS OUTPUT SelectionSummary=_removed ParameterEstimates=_estimate;
      PROC GLMSELECT DATA = &dsn;
         class &cvar;
         model &outcome = &var/selection=backward include=&inc slstay=&slstay select=SL 
          MAXSTEP=1 hierarchy=single;
         %if &weight ~= %STR() %then %do; 
            weight &weight;
         %end;
      RUN;

      /* PROC GLMSELECT creates the removed data set whether or not any variables are removed */
      /* If no variables were removed then the selection process is done */
      /* Get name of variable removed */
      PROC SQL noprint;
         select UPCASE(EffectRemoved), count(*)
         into :remove, :continue
         from _removed
         /* Exclude full model step */
         where EffectRemoved ~= ' ';
      quit;

      /* Update variable list with selected vars only */
      /* Note that the order of variables in an interaction term can be revsered so the */
      /* method used to update the categorical variable list cannot be used.  The method */
      /* below overcomes this problem. */

      /* Save order */
      DATA _est2;
         set _estimate;
         order = _n_;
      RUN;

      /* Get unique list of variable names */
      PROC SORT DATA = _est2 nodupkey;
         by Effect;
         where Effect ~= 'Intercept';
      RUN;

      /* Return to original order */
      PROC SORT DATA = _est2;
         by order;
      RUN;

      /* Update variable list */
      PROC SQL noprint;
         select UPCASE(effect) into: var separated by ' '
         from _est2;
      QUIT;

      %if &continue = 1 %then %do;

         %let remlist = &remlist &remove;

         %if &debug = T %then %do;
            %put remlist &remlist;
         %end;

         /* Update class variable list */
         %let newcvar =;
         %do i = 1 %to &cvar_cnt;

            %let c&i = %SCAN(&cvar, &i, ' ');              
            /* Remove from list if it is the removed variable */
            /* Rescan to just get variable name not order as well */
            %if %BQUOTE(%SCAN(&&c&i, 1, %STR( ))) ~= %BQUOTE(&remove) %then %do;
               /* Recombine individual variables into a new list */
               %let newcvar = &newcvar &&c&i; 
            %end;
         %end;
 
         /* Reset var list */
         %let cvar = &newcvar;
         /* Trim extra blanks */
         %let cvar = %SYSFUNC(TRANWRD(&cvar,%STR(  ),%STR()));

         /* Update variable counts */
         %let var_cnt = %EVAL(&var_cnt-1);
         /* Count number of class variables */
         %let cvar_cnt = %sysfunc(countw(&cvar,' ')); 

         /* Delete dataset */
         proc datasets lib=work;
            DELETE _removed;
         quit;  
      %end;

   %end;

   /* Save selected vars in macro variables */
   %let _finalvar = &var;
   %let _finalcvar = &cvar;

   /* Produce report of final model */
   %if &report = T %then %do;

      /* If any variables were removed from the model create a list */
      %if &remlist ~= %STR() %then %do;

         /* Get variable labels for footnote */
         PROC CONTENTS DATA = &dsn (keep=&remlist) out=_cont;
         RUN;

         DATA _cont;
            set _cont end=last;
             if label = ' ' then label = name;
             if _N_ ~= 1 and last then label = 'and ' || label;
         RUN;

         PROC SQL noprint;
            select label into :remLab separated by ', '
             from _cont;
         QUIT;

         %let foottext = The following variables were removed from the model: &remLab..;

      %end;
      %else %let foottext = No variables were removed from the model.;

      /* Remove outer quotes before macro call */
      %if %sysevalf(%superq(event)~=,boolean) %then %do;
         %let event = %sysfunc(DEQUOTE(&event));
      %end;

      /* Fit final model */
      ODS OUTPUT NObs = nobs 'Type III Model ANOVA'=type3 ParameterEstimates=estimate;
      PROC GLM DATA = &dsn ORDER=&order;
         class &_finalcvar;
         model &outcome = &_finalvar/solution clparm;
         %if &weight ~= %STR() %then %do; 
            weight &weight;
         %end;
      RUN;
      QUIT;

      /* Produce report */
      %MULTIPLE_LINREG(DSN=&dsn,
          OUTPATH=&outpath,
          FNAME=&filename,
          FOOTNOTE="Backward selection with an alpha level of removal of &slstay was used.  &foottext",
          TYPE3 = &type3,
          debug=&debug);

   %end;

   /* Only delete files if not in debug mode */
   %if &debug ~= T %then %do;
      proc datasets lib=work memtype=data noprint;  
         delete %if &report = T %then %do; estimate nobs type3 %end; _removed _estimate 
          %if &report = F %then %do; _est2 %end;;
      quit;  
   %end;

   /* Final variables selected (get rid of double spaces */
   %put Categorical variables selected: &_finalcvar;
   %put All variables selected: &_finalvar;

%mend linreg_sel;

/*********************************************************************************************/
