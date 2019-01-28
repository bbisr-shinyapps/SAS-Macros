/**********************************************************************************************
***********************************************************************************************
Module Name: COXPHFIRTH
Created Date/Author/Contact: June 10, 2014/Dana Nickleach/dnickle@emory.edu
Last Update Date/Person/Contact: 
Current Version: V1

Purpose:  To conduct univariate Cox proportional hazards regression using Firth's bias 
correction and penalized likelihood confidence intervals and p-values by calling R.  

Notes:  R must be installed on your computer as well as the package coxphf.  This module will
be used by UNI_PHREG V23.  The data sets _firth and _firthnames will be produced containing
the results.  _firth will contain the estimates and _firthnames will contain the names of the
predictors and levels in the case of a categorical variable.

Parameters: 
DSN            The name of the data set to be analyzed.
EVENT          The variable name for time to event outcome.
CENSOR         The variable name for censoring indicator. It is required that 1 is used for the 
               event and 0 for censored cases. 
VAR            The name of the categorical or numeric covariate. 
REF            The reference category of the categorical variable.  This should be set to "NA"
               for numeric covariates.
 
**********************************************************************************************/
/* COXPHF Module */
PROC IML;

   start coxphfirth(dsn,event,censor,var,ref);

   /* Export Data to R*/

   RUN ExportDataSetToR(dsn, "data");

   SUBMIT e=event c=censor x=var r=ref/ R;

      if ("&r"!="NA"){
         #Set reference level for categorical variables
         data$&x <- relevel(data$&x, ref="&r")
      }

      require(coxphf)
      #Remove missing values
      datanew = data[complete.cases(data$&x,data$&e,data$&c),]
      result=coxphf(datanew, formula=Surv(&e,&c)~&x)
      #Get HR, CI, and p-value
      report=cbind(HRest=exp(result$coefficient),LCL=result$ci.lower,UCL=result$ci.upper,pval=result$prob,nobs=result$n)
      #Variable names 
      names=row.names(report)
   ENDSUBMIT;

   /* Transfer confidence intervals from R to SAS */
   run ImportMatrixFromR(firth, "report" );
   run ImportMatrixFromR(names, "names" );

   /* Since firth is numeric and names is character they can't be combined into one matrix in IML */
   varNames = {"HR", "LCL", "UCL", "pval", "nobs"};

   /* Create a SAS data set */
   create _firth from firth [colname=varNames]; 
   append from firth;

   /* Create a SAS data set */
   create _firthnames from names [colname="var"]; 
   append from names;

   finish;
   store module=coxphfirth;
QUIT;
