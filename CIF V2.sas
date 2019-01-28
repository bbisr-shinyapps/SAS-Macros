 /*****************************************************************/
 /*          S A S   A U T O C A L L   L I B R A R Y              */
 /*                                                               */
 /*   NAME : CIF                                                  */
 /*   TITLE: Estimate Cumulative Incidence Function and           */
 /*          Compare Cumulative Incidence Functions across Samples*/
 /* PRODUCT: STAT                                                 */
 /*  SYSTEM: ALL                                                  */
 /*    KEYS: survival data, failure probability, competing risks  */
 /*   PROCS:                                                      */
 /*    DATA:                                                      */
 /* SUPPORT: Guixian Lin                UPDATE:  Mar2012          */
 /*     REF: Klein  and Moeschberger   (2003)                     */
 /*          Marubini and Valsecchi    (1995)                     */
 /*          Gray                      (1988)                     */
 /*    MISC: The %CIF macro requires Base SAS and SAS/IML.        */
 /*          SAS/GRAPH is also required in SAS 9.2                */
 /*****************************************************************/
/**********************************************************************************************

Revision History:

Version		Date		By		Description
V2			9/20/12		DCN		Added Gray's p-value to the plot, added yaxis statement, and 
								added NOPRINTTEST to options.
V2			11/29/12	DCN		Removed nonotes option.

**********************************************************************************************/

/*------------------------------------------------------------------
The %CIF macro computes the crude cumulative incidence function
estimates for homogeneous (no covariates) survival data whose
endpoints are subjected to competing risks (Klein and
Moeschberger, 2003). Standard errors and pointwise confidence
limits are also computed. The estimated cumulative incidence curve
is plotted as a  step function using ODS graphics. You can also
display the pointwise confidence limits in the plot by specifying
the PLOTCL keyword in the OPTIONS= argument. If you do not want any
plots, you can suppress it by specifying the NOPLOT keyword in the
OPTIONS= argument.

If the data contains multiple samples identified by the GROUP= variable,
the %CIF macro computes a cumulative incidence curve for each group, and
the macro also compares cumulative incidence functions between groups
using the method of Gray (1988).


--------------------------------------
EXAMPLES:

-----Example 1:

Consider the data in Table 10.1 of Marubini and
Valsecchi (1995). The data include two treatment groups (A and B), and
two competing events (1 and 2). The variable status has 3 distinct
values: value 0 indicates a censored time, and values 1 and 2 represent
two competing events.


Here we illustrate how to estimate the cumulative incidence function within
each group with the event of interest(=2), and test whether there is any
difference  between these two cumulative incidence functions.


data a (drop=_i_);
   array type1(10) _temporary_ (1 13 17 30 34 41 78 100 119 169);
   array type2(20) _temporary_ (1 6 8 13 13 15 33 37 44 45
                                63 80 89 89 91 132 144 171 183 240);
   array type3(5)  _temporary_ (34 60 63 149 207);
   do _i_=1 to 10;
      Status=1;
      t=type1{_i_};
      output;
      end;
   do _i_=1 to 20;
      Status=2;
      t=type2{_i_};
      output;
      end;
   do _i_=1 to 5;
      Status=0;
      t=type3{_i_};
      output;
      end;
run;

data b (drop=_i_);
   array type1(10) _temporary_ (7 16 16 20 39 49 56 73 93 113);
   array type2(20) _temporary_ (1 2 4 6 8 9 10 13 17 17
                                17 18 18 27 29 39 50 69 76 110);
   array type3(5)  _temporary_ (34 60 63 78 149);
   do _i_=1 to 10;
      Status=1;
      t=type1{_i_};
      output;
      end;
   do _i_=1 to 20;
      Status=2;
      t=type2{_i_};
      output;
      end;
   do _i_=1 to 5;
      Status=0;
      t=type3{_i_};
      output;
      end;
run;

data ab;
   label t="Weeks";
   set a(in=ina) b;
   if ina then group='A';
   else        group='B';

run;



%CIF(data=ab, out=est, time=t, status=status, group=group, event=2, title=CIF curves, title2=for Status=2);




-----Example 2:
Consider the data on 137 Bone Marrow Transplant Patients in Appendix D of Klein
and Moeschberger (2003). In this data, Patients are grouped into 3 categories:
ALL (38 patients), AML low risk (54 patients), and AML high risk (45). The competing
events are relapse and death in remission.  The variables used in the
example are:
Gender: patient gender (1-Male, 0-Female)
Status: event indicator (0=censored, 1=death in remission, 2=relapse)
   Grp: disease categories (1-ALL, 2-AML low risk, 3-AML high risk)
 Ftime: failure time.

The following DATA step creates the SAS data set BMT. Here we use the data set BMT
to illustrate how to compute the cumulative incidence functions with the event of interest
for different groups, and test whether there is any difference between the CIF across
groups using stratified data by specifying the strata= variables.


proc format;
   value grpLabel 1='ALL' 2='AML low risk' 3='AML high risk';
run;

data BMT;
        input Grp Ftime Status Gender@@;
        label Ftime="Days";
        format Diagnosis grpLabel.;
datalines;
1       2081       0       1       1       1602    0       1
1       1496       0       1       1       1462    0       0
1       1433       0       1       1       1377    0       1
1       1330       0       1       1       996     0       1
1       226        0       0       1       1199    0       1
1       1111       0       1       1       530     0       1
1       1182       0       0       1       1167    0       0
1       418        2       1       1       383     1       1
1       276        2       0       1       104     1       1
1       609        1       1       1       172     2       0
1       487        2       1       1       662     1       1
1       194        2       0       1       230     1       0
1       526        2       1       1       122     2       1
1       129        1       0       1       74      1       1
1       122        1       0       1       86      2       1
1       466        2       1       1       192     1       1
1       109        1       1       1       55      1       0
1       1          2       1       1       107     2       1
1       110        1       0       1       332     2       1
2       2569       0       1       2       2506    0       1
2       2409       0       1       2       2218    0       1
2       1857       0       0       2       1829    0       1
2       1562       0       1       2       1470    0       1
2       1363       0       1       2       1030    0       0
2       860        0       0       2       1258    0       0
2       2246       0       0       2       1870    0       0
2       1799       0       1       2       1709    0       0
2       1674       0       1       2       1568    0       1
2       1527       0       0       2       1324    0       1
2       957        0       1       2       932     0       0
2       847        0       1       2       848     0       1
2       1850       0       0       2       1843    0       0
2       1535       0       0       2       1447    0       0
2       1384       0       0       2       414     2       1
2       2204       2       0       2       1063    2       1
2       481        2       1       2       105     2       1
2       641        2       1       2       390     2       1
2       288        2       1       2       421     1       1
2       79         2       0       2       748     1       1
2       486        1       0       2       48      2       0
2       272        1       0       2       1074    2       1
2       381        1       0       2       10      2       1
2       53         2       0       2       80      2       0
2       35         2       0       2       248     1       1
2       704        2       0       2       211     1       1
2       219        1       1       2       606     1       1
3       2640       0       1       3       2430    0       1
3       2252       0       1       3       2140    0       1
3       2133       0       0       3       1238    0       1
3       1631       0       1       3       2024    0       0
3       1345       0       1       3       1136    0       1
3       845        0       0       3       422     1       0
3       162        2       1       3       84      1       0
3       100        1       1       3       2       2       1
3       47         1       1       3       242     1       1
3       456        1       1       3       268     1       0
3       318        2       0       3       32      1       1
3       467        1       0       3       47      1       1
3       390        1       1       3       183     2       0
3       105        2       1       3       115     1       0
3       164        2       0       3       93      1       0
3       120        1       0       3       80      2       1
3       677        2       1       3       64      1       0
3       168        2       0       3       74      2       0
3       16         2       0       3       157     1       0
3       625        1       0       3       48      1       0
3       273        1       1       3       63      2       1
3       76         1       1       3       113     1       0
3       363        2       1
;
run;


%CIF(data=bmt,time=Ftime, status=Status, group=Diagnosis, event=1 2, censored=0, strata=Gender, alpha=.05);



----------------------------------------------------------------------
References:
1. Marubini, E. and Valsecchi, M.G. (1995), Analysing survival data from clinical trials
   and observational studies, John Wiley.
2. Gray, R.J. (1988). "A class of K-sample tests for comparing the cumulative incidence of
   a competing risk," Annals of statistics, 16(3), 1141--1154.
3. Klein, J.P. and Moeschberger, M.L., (2003), Survival analysis: techniques for censored
   and truncated data, Springer Verlag.





----------------------------------------------------------------------

 DISCLAIMER:
 THIS INFORMATION IS PROVIDED BY SAS INSTITUTE INC. AS A SERVICE
 TO ITS USERS.  IT IS PROVIDED "AS IS".  THERE ARE NO WARRANTIES,
 EXPRESSED OR IMPLIED, AS TO MERCHANTABILITY OR FITNESS FOR A
 PARTICULAR PURPOSE REGARDING THE ACCURACY OF THE MATERIALS OR CODE
 CONTAINED HEREIN.

----------------------------------------------------------------------*/


%macro CIF(
   data=,           /* Name of the input data set                             */
   out=cifEstimate, /* Name of the output data set that                       */
                    /* contains the cumulative incidence estimates            */
   time=Time,       /* Time variable                                          */
   status=Status,   /* Numeric variable whose values indicate  whether an     */
                    /* observation corresponds to the event of interest, the  */
                    /* competing event, or is censored                        */
   event=1,         /* Values of status= variable that correspond to the event*/
                    /* interest. Multiple values are separated by blanks      */
   censored=0,      /* Values of status= variable that correspond to censored */
                    /* observations. Multiple values are separated blanks     */
   group=,          /* Variable that defines the groups for comparison        */
   strata=,         /* Variable that defines the  strata for stratified test  */
   alpha=.05,       /* Complement of the confidence level for pointwise       */
                    /*confidence limits                                       */
   rho=0,           /* Power of the weight function in the test               */
   SE=1,            /* Method to compute the standard error of cumulative     */
                    /* incidence functions. 1=delta method;2=counting process */
   TITLE=,          /* First title for the cumulative incidence function      */
                    /* plot (without quote)                                   */
   TITLE2=,         /* Second title for the cumulative incidence function     */
                    /* plot (without quote)                                   */
   options=,        /* Multiple keywords are separated by blanks              */
                    /* NOPLOT suppresses plotting the cumulative incidence    */
                    /*        curves                                          */
                    /* PLOTCL plots the pointwise confidence limits for       */
                    /*        cumulative incidence curves                     */
                    /* NOCIFEST suppresses the CIF estimate table             */
                    /* NOTEST suppresses the test for differences among       */
                    /*        cumulative incidence functions                  */
   				    /* NOPRINTTEST suppresses the printing of the table       */
					/* 		  containing the results from the test.           */
                    /*                                                        */
   );;


%local Confidence savenote noplotflag plotclflag notestflag nocifflag nev s11 inconsistentgrpnumflag;
%let savenote     = %sysfunc(getoption(notes));

/*options nonotes;*/


/* Check if time= and status= variables okey */
proc contents data=&data noprint out=_c; run;

%let TimeVar=1;
%let StatusVar=1;
%let GroupVar=1;
%let StrataVar=1;
%let nocifflag=0;
%let inconsistentgrpnumflag=0;

data _null_;
   length x y $32;
   set _c;
   x=lowcase(NAME);
   y=lowcase("&time");
   if x=y then do;
      if TYPE=1 then call symput('TimeVar', put(0,1.));
      else           call symput('TimeVar', put(2,1.));
   end;
   y=lowcase("&status");
   if x=y then do;
      if TYPE=1 then call symput('StatusVar', put(0,1.));
      else           call symput('StatusVar', put(2,1.));
   end;
   %if ((&group^=) ) %then %do;
           y=lowcase("&group");
           if x=y then do;
                  if TYPE=1 then call symput('GroupVar', put(0,1.));
                  else           call symput('GroupVar', put(2,1.));
                  end;
           %end;
   %else
       call symput('GroupVar', put(3,1.));;

   %if ((&strata^=) ) %then %do;
           y=lowcase("&strata");
           if x=y then do;
                  if TYPE=1 then call symput('StrataVar', put(0,1.));
                  else           call symput('StrataVar', put(2,1.));
           end;
           %end;
   %else
                call symput('StrataVar', put(3,1.));;

run;


%if &TimeVar=1 %then %do;
   %let _xrc_=
     %str(Variable &time not found in the input data set);
   %put ERROR: &_xrc_..;
   %goto exit;
%end;
%if &TimeVar=2 %then %do;
   %let _xrc_=
     %str(Variable &time is not a numeric variable);
   %put ERROR: &_xrc_..;
   %goto exit;
%end;
%if &StatusVar=1 %then %do;
   %let _xrc_=
    %str(Variable &status not found in the input data set);
   %put ERROR: &_xrc_..;
   %goto exit;
%end;
%if &StatusVar=2 %then %do;
   %let _xrc_=
    %str(Variable &status is not a numeric variable);
   %put ERROR: &_xrc_..;
   %goto exit;
%end;

%if &GroupVar=1 %then %do;
   %let _xrc_=
    %str(Variable &group not found in the input data set);
   %put ERROR: &_xrc_..;
   %goto exit;
%end;

%if &StrataVar=1 %then %do;
   %let _xrc_=
    %str(Variable &strata not found in the input data set);
    %put ERROR: &_xrc_..;
   %goto exit;
%end;



data _null_;
   /*----- check OPTIONS= -----*/
   list= lowcase(symget('options'));
   i = index(list, 'noplot'); if i then substr(list, i, 6) = ' ';
   call symput('noplotflag', put(i gt 0, 1.));
   i = index(list, 'plotcl'); if i then substr(list, i, 6) = ' ';
   call symput('plotclflag', put(i gt 0, 1.));
   i = index(list, 'notest'); if i then substr(list, i, 6) = ' ';
   call symput('notestflag', put(i gt 0, 1.));
   i = index(list, 'nocifest'); if i then substr(list, i, 8) = ' ';
   call symput('nocifestflag', put(i gt 0, 1.));
   i = index(list, 'noprinttest'); if i then substr(list, i, 11) = ' ';
   call symput('noprinttestflag', put(i gt 0, 1.));

   /*----- Confidence Level -------*/
   alpha=.05;
   %if &alpha ne %then %do; alpha=&alpha; %end;
   call symput('Confidence',trim(put(100*(1-alpha),best8. -l)));
   /*---- No of event= and compete= values ----*/
   ne=countw(symget('event'),' ');
   call symput('ne',put(ne,best8. -l));
   ncens=countw(symget('censored'),' ');
   call symput('ncens',put(ncens,best8. -l));

run;


%if &data= %then %do;
   %let _xrc_= %str(The DATA= parameter is not specified);
   %put ERROR: &_xrc_..;
   %goto exit;
   %end;
%else
   %do;

   data _null_;
   set &data;
   array event[&ne] _temporary_(&event);
   array censored[&ncens] _temporary_(&censored);;
   zevent=0;
   zcensored=0;
   do i=1 to &ne;
     if &status = event[i] then do;
       zevent=1;
       goto skip1;
     end;
   end;
   skip1:;

   do i=1 to &ncens;
      if &status = censored[i] then do;
         zcensored=1;
         go to skip2;
      end;
   end;
   skip2:;




   if zevent=1 and zcensored=1 then do;
      call symput('s11',trim(left(put(&status,best8.))));
      stop;
   end;
   if zevent=1 then do;
       call symput('nev',put(1,1.));

   end;


   run;


%if &nev eq %then %do;
  %let _xrc_=
    %str(The total number of EVENT is 0);
  %put ERROR: &_xrc_..;
  %goto exit;
%end;

%if &s11 ne %then %do;
  %let _xrc_=
    %str(&Status=%sysfunc(trim(&s11)) is both an EVENT= and CENSORED= value);
  %put ERROR: &_xrc_..;
  %goto exit;
%end;
run;

data _null_;
   set &data;
   length l $ 60;
   call vname(&time, l);
   call symputx('tname', l, 'L');
   call label(&time, l);
   call symputx('tlabel', l, 'L');
   %if ((&group^=) ) %then %do;
      call vname(&group, l);
      call symputx('gname', l, 'L');
      call label(&group, l);
      call symputx('glabel',l , 'L');
      call symputx('grpFmt', VFORMAT(&group));
      %end;
   %else %do;
      call symputx('glabel',"Group" , 'L');
          call symputx('gname', "Group", 'L');;
      %end;

    %if ((&strata^=) ) %then %do;
      call vname(&strata, l);
      call symputx('stratname', l, 'L');
      call label(&strata, l);
      call symputx('stratlabel',l , 'L');
      call symputx('stratFmt', VFORMAT(&strata));
      %end;
   %else %do;
      call symputx('stratlabel',"Strata" , 'L');
          call symputx('stratname', "Strata", 'L');;
      %end;


run;

/*create a data set by ignoring those missing observations*/
data &data.1(drop=&group);
   set &data;
   %if ((&group^=) ) %then %do;
      f&group= put(&group, &grpFmt);
      IF CMISS(&time,&status, f&group) > 0 THEN DELETE;
      %end;
   %if ((&strata^=) ) %then %do;
      f&strata= put(&strata, &stratFmt);
      IF CMISS(&time, &status, f&strata) > 0 THEN  DELETE;
      %end;
   %end;


run;


data _null_;
   call symputx('nobs0',nobs0);
   stop;
   set &data nobs=nobs0;
run;

data _null_;
   call symputx('nobs1',nobs1);
   stop;
   set &data.1 nobs=nobs1;
run;


%if &nobs0>&nobs1 %then %do;
   %put NOTE: There are %eval(&nobs0-&nobs1) missing observations;

%end;


proc iml;
   use &data.1;
   read all var{&time} into time;
   read all var{&status} into status;
   nobs= nrow(time);
   %if %bquote(&group)= %then %do;
     group= j(nobs,1,1);
   %end;
   %else %do;
   read all var{f&group} into group;
   %end;


   %if %bquote(&strata)= %then %do;
     strat= j(nobs,1,1);
   %end;
   %else %do;
     read all var{f&strata} into strat;
   %end;


   close &data.1;

* subroutine to calculate the cumulative incidence and variance;
start estimate_cif(y, causecode) global(time, cif, vcif1, vcif2, n2);
/*************************************************************
* y is the failure times (sorted in increasing order);
* ic is the indicator for censoring
*   (= 0 if censored, =1 if failed from any cause);
* icc is the indicator for failure from the cause of interest
*   (= 0 if failed from cause of interest, =1 otherwise);
* n is the # of observations (length of y);
* n2 = n+2:  the number of the unique time points from cause 1 with 2 extra points (0 and largest time point within this group);
* ut is the unique failure times from the cause of interest;
* cif is the estimated cumulative incidence function;
* vcif1 and vcif2 are the estimated variance of cif from Marubini method and Gray method respectively;
**************************************************************/
   surv_km= 1;
   v11= 0; v21= 0;
   v12= 0; v22= 0;
   v13= 0; v23= 0;


   cif_pos= 1;
   n= nrow(y);
   nrisk=n;

   ut=unique(y)`;
   n_ut=nrow(ut);

   cif= j(n2,1,0);
   vcif1= j(n2,1,0);
   vcif2= j(n2,1,0);


   do i=1 to n_ut;
      *compute the number of events of interest and competing events;
      nevent=sum(causecode[loc(y=ut[i])]=1);
      ncompetevent=sum(causecode[loc(y=ut[i])]=2);
      n_allevents = nevent+ ncompetevent;

      *compute survival function and cumulative incidence function;
      if n_allevents^= 0 then do; *overall survival changes whenever a failure occurs;
         if (i>1)  then nrisk =n - min(loc(y=ut[i]))+ 1 ;
         surv_km_n= surv_km* (nrisk- n_allevents)/ nrisk;
         if nevent> 0 then do;
            cif_pos= cif_pos + 1;
            cif[cif_pos]= cif[cif_pos- 1]+ surv_km * nevent/ nrisk;
            end;


            * variance estimator using delta method- (10.12) in Marubini and Valsecchi (1995);
            * Analysing Survival Data from Clinical Trials and Observational Studies;
            t1=0;
            if (nrisk- n_allevents>0) then         t1= n_allevents/ nrisk/ (nrisk- n_allevents);
            t2= surv_km* nevent/ nrisk/ nrisk;
            t3= surv_km* (1- nevent/ nrisk)* t2;
            t4= cif[cif_pos]* cif[cif_pos];
            v11= v11+ t1;
            v12= v12- 2* (cif[cif_pos]* t1+ t2);
            v13= v13+ t4* t1+ 2* cif[cif_pos]* t2 + t3;
            vcif1[cif_pos]= t4* v11+ cif[cif_pos]* v12+ v13;



            *variance estimator using counting process- Gray (1988);

            if (ncompetevent> 0 & surv_km_n> 0) then do;
               vt1=1;
               if( ncompetevent > 1) then vt1= 1- (ncompetevent - 1)/(nrisk- 1);
               vt2= surv_km* surv_km* vt1* ncompetevent/(nrisk* nrisk);
               vt3= 1/ surv_km_n;
               vt4= cif[cif_pos]/ surv_km_n;
               v21= v21+ vt2* vt4* vt4;
               v22= v22+ vt2* vt3* vt4;
               v23= v23+ vt2* vt3* vt3;
               end;
            if (nevent> 0) then do;
               vt1=1;
               if( nevent> 1 ) then vt1= 1-(nevent- 1)/(nrisk- 1);
               vt2= surv_km* surv_km* vt1* nevent/(nrisk* nrisk);
               if surv_km_n> 0 then vt3= 1/ surv_km_n; else vt3= 0;
               vt4= 1+ vt3* cif[cif_pos];
               v21= v21+ vt2* vt4* vt4;
               v22= v22+ vt2* vt3* vt4;
               v23= v23+ vt2* vt3* vt3;
               vcif2[cif_pos]= v21 - 2*cif[cif_pos]*v22 +  cif[cif_pos]*cif[cif_pos]*v23;
               end;

            surv_km= surv_km_n;

         end;


   end;       /*end do i=*/


   if cif_pos< n2 then do;
      /*deal with the largest time point with the group*/
      cif_pos= cif_pos+ 1;
      time[cif_pos]= ut[n_ut];
      cif[cif_pos]= cif[cif_pos- 1];
      vcif1[cif_pos]= vcif1[cif_pos- 1];
      vcif2[cif_pos]= vcif2[cif_pos- 1];

   end;
   if (sum(vcif1<0)>0) then   vcif1[loc(vcif1<0)]=0;
   if (sum(vcif2<0)>0) then   vcif2[loc(vcif2<0)]=0;


   vcif1= sqrt(vcif1);
   vcif2= sqrt(vcif2);

finish estimate_cif;



/* subroutine for calculating k-sample test statitsics and variance within stratum;
* for competing risk for Gray's method (1988);*/
start test_ksample(y, causecode, ig, rho) global(score, vcov);
*  y is the failure times (sorted in ascending order)
*  causecode is cause indicator
*   (= 0 if censored, =1 if failed from the cause of interest,
*       2 if failed from some other causes.)
*  ig denotes group membership

*  rho is the power of the weight function in the test statistic


   *  ng denotes the number of groups, ng1= ng- 1;
   ng= ncol(unique(ig));
   * risk set size in each at the current failure time (= sample size at 0);
   nrisk= j(ng,1,0);

   do i= 1 to ng;
      nrisk[i]= nrisk[i] + sum(ig=i);
   end;

   ng1=ng-1;
   ng2= ng* ng1/ 2;
   score= j(ng1,1,0);
   vcov= j(ng2,1,0);

   l= 0;
   cifgleft= j(ng,1,0);
   cifg=  j(ng,1,0);
   cif0left= 0;
   cif0= 0;

   vtvec=  j(ng,1,0);
   skmgleft= j(ng,1,1);
   skmg= j(ng,1,1);


   vtmatrix= j(ng1,ng,0);
   c= j(ng,ng,0);


   ll= 1;

   ut=unique(y)`;
   n_ut=nrow(ut);

   do i=1 to n_ut;
         d= j(3,ng,0);
         yloc=loc(y=ut[i]);
      do ii= yloc[1] to yloc[ncol(yloc)];
         k= causecode[ii]+ 1;   * causes code:0, 1, 2;
         d[k,ig[ii]]= d[k,ig[ii]] + 1;
      end;
      nd1= d[2,+];
      nd2= d[3,+];

      if (nd1^= 0 | nd2^= 0) then do;
         *hdot is sum(hhat), Rdot is sum(R), skmgleft is Shat, nrisk is Y, cifgleft is CIF;
         hdot= 0; Rdot= 0;
         do g= 1 to ng;
            if (nrisk[g]> 0) then do;
               td= d[2,g] + d[3,g];
               skmg[g]= skmgleft[g]* (nrisk[g]- td)/ nrisk[g];
               cifg[g]= cifgleft[g]+ skmgleft[g]* d[2,g]/ nrisk[g];
               hdot= hdot+ nrisk[g]/ skmgleft[g];
               Rdot= Rdot+ nrisk[g]* (1 - cifgleft[g])/ skmgleft[g];
               end;
         end;

         cif0= cif0left+ nd1/hdot;
         grho= (1- cif0left)**rho;

         * score statistic;
         do g=1 to ng1;
            if (nrisk[g]> 0) then score[g]= score[g]+ grho* (d[2,g]- nd1* nrisk[g]* (1- cifgleft[g])/ skmgleft[g]/ Rdot);
         end;

         a= j(ng, ng,0);
         do g= 1 to ng;
            if (nrisk[g]> 0) then do;
               t1= nrisk[g]/ skmgleft[g];
               a[g,g]= grho* t1* (1- t1/ hdot);
               c[g,g]= c[g,g]+ a[g,g]* nd1/hdot/(1- cif0left);
               k= g + 1;
               *if (k<= ng) then do j= k to ng;
               do j= g+1 to ng;
                  if (nrisk[j]> 0) then do;
                     a[g,j]= - grho* t1* nrisk[j]/ skmgleft[j]/ hdot; a[j,g]= a[g,j];
                     c[g,j]= c[g,j]+ a[g,j]* nd1/ hdot/ (1 - cif0left); c[j,g]= c[g,j];
                     end;
               end;
               end;
         end;


         * variance estimators;
         if nd1> 0 then do k=1 to ng;
            if (nrisk[k]> 0) then do;
               if skmg[k]> 0 then vt1= 1- (1- cif0)/ skmg[k]; else vt1= 1;
               if nd1> 1 then vt2= 1- (nd1- 1)/ (hdot* skmgleft[k]- 1); else vt2= 1;
               vt3= vt2* skmgleft[k]* nd1/ (hdot* nrisk[k]);
               vtvec[k]= vtvec[k]+ vt1* vt1* vt3;
               do g= 1 to ng1;
                  vt4= a[g,k]- vt1* c[g,k];
                  vtmatrix[g,k]= vtmatrix[g,k]+ vt4* vt1* vt3;
                  do j=1 to g;
                     l= g* (g- 1)/ 2+ j;
                     vt5= a[j,k]- vt1* c[j,k];
                     vcov[l]= vcov[l]+ vt3* vt4* vt5;
                  end;
               end;
               end;
            end;

         if (nd2^= 0) then do k=1 to ng;
            if (skmg[k]> 0 & d[3,k]> 0) then do;
               vt1= (1- cif0)/ skmg[k];
               vt2= 1;
               if (d[3,k]> 1) then vt2= 1- (d[3,k] - 1)/ (nrisk[k]- 1);
               vt3= vt2* ((skmgleft[k]** 2)* d[3,k])/ (nrisk[k]** 2);
               vtvec[k]= vtvec[k]+ vt1* vt1* vt3;
               do g=1 to ng1;
                  vt4= vt1* c[g,k];
                  vtmatrix[g,k]= vtmatrix[g,k]- vt4* vt1* vt3;
                  do j=1 to g;
                     l= g* (g- 1)/ 2+ j;
                     vt5= vt1* c[j,k];
                     vcov[l]= vcov[l]+ vt3* vt4* vt5;
                  end;
               end;
               end;
            end; /*if (nd2^= 0)*/

         end; /*if (nd1^= 0 | nd2^= 0)*/

      *update the risk sets and the index of the next failure time;

      do ii= yloc[1] to yloc[ncol(yloc)];
         nrisk[ig[ii]]= nrisk[ig[ii]] - 1;
      end;
      cif0left= cif0;
      cifgleft= cifg;
      skmgleft= skmg;



   end;  /* end do i=*/

   pos= 0;
   do g= 1 to ng1;
      do j= 1 to g;
         pos= pos+ 1;
         do k= 1 to ng;
            vcov[pos]= vcov[pos]+ c[g,k]* c[j,k]* vtvec[k];
            vcov[pos]= vcov[pos]+ c[g,k]* vtmatrix[j,k];
            vcov[pos]= vcov[pos]+ c[j,k]* vtmatrix[g,k];
         end;
      end;
   end;

finish test_ksample;





   * subroutines end and the calculation starts here;
   * status need to be numerical value;

   x= time|| status;

   *call sort(x, {1 2});
   call sortndx(sidx, x, {1 2});
   x=x[sidx,];



   group=group[sidx];
   strat=strat[sidx];
   ugg= unique(group);
   ng= ncol(ugg);
   call symput('ng',left(char(ng)));

   * censind= 0 if censored, =1 if failed;

   censList={&censored}`;
   eventList={&event}`;
   censind=j(nrow(x),1,0);

   eventind=censind;
   do i= 1 to &ncens;
     censind= censind |(x[,2]= censList[i]);
   end;
   do i= 1 to &ne;
     eventind= eventind |(x[,2]= eventList[i]);
   end;

   causecode= 2* (^censind)- eventind;


   * prepare the data;
   x= x|| censind|| causecode;
   if ng> 1 then do;
      ng1= ng - 1;
      ng2= ng* ng1/ 2;
      vv= j(ng1,ng1,0); *the final covariance matrix;
      vcov_st= j(ng2,1,0);  vcov= vcov_st; *vcov_st for each stratum, vs add v across strata ;
      score_st= j(ng1,1,0); score= score_st; *s is test statiatics by stratum,ss add s across strata;

   end;

   ust= unique(strat);
   nstrata= ncol(ust);
   call symput('nstrata',left(char(nstrata)));

  /* calculate CIF and  variance*/
   if (&nocifflag=0) then do;

          do k= 1 to nstrata;

                  do g= 1 to ng;
                         gind= (group= ugg[g]);
                         *ncg= sum(gind);
                         * subset of data for group g and stratum k;
                         grpstra=(gind= 1) & (strat= ust[k]) & (x[,4]=1);

                         if (sum(grpstra)>0) then do;

                                 xsub= x[loc((gind= 1) & (strat= ust[k])), ];

                                 * aa are failure times for cause 1, n2 is the length of CIF;
                                 aa= xsub[loc(xsub[,4]= 1), 1];
                                 n2 = ncol(unique(aa))+ 2;

                                 run estimate_cif(xsub[,1], xsub[,4]);

                                 time1=0//(unique(aa))`//xsub[nrow(xsub),1];
                                 tmp= time1|| cif|| vcif1|| vcif2;  *time,CI,variance,group;
                                 if (g= 1 & k=1) then do;
                                        cifEST= tmp;
                                        grpLabel=j(n2,1,ugg[g]);
                                        if(nstrata>1) then stratLabel=j(n2,1,ust[k]);

                                        end;
                                 else  do;
                                        cifEST= cifEST//tmp;
                                        grpLabel= grpLabel//j(n2,1,ugg[g]);
                                        if(nstrata>1) then stratLabel=stratLabel//j(n2,1,ust[k]);
                                        end;
                                end; /*(sum(grpstra)>0)*/
                  end;
           end;



      /*create SAS data set for CIF estimate*/
      cifloc=loc(cifEST[,2]>0);
      ncif= nrow(cifEST);
      lowci=j(ncif, 1, 0);
      upci =j(ncif, 1, 0);
      z=probit(1-&alpha/2.0);

      /*confidence intervals:  log-log tranform*/
      lowci[cifloc] = cifEST[cifloc, 2]##exp((-z*cifEST[cifloc, 3]) /(cifEST[cifloc,2]#LOG(cifEST[cifloc,2])));
      upci[cifloc] = cifEST[cifloc, 2]##exp((z*cifEST[cifloc, 3]) /(cifEST[cifloc,2]#LOG(cifEST[cifloc,2])));

      if (&SE=1) then do;
         /*Counting process method -- Gray*/
         lowci[cifloc] = cifEST[cifloc, 2]##exp((-z*cifEST[cifloc, 4]) /(cifEST[cifloc,2]#LOG(cifEST[cifloc,2])));
         upci[cifloc] = cifEST[cifloc, 2]##exp((z*cifEST[cifloc, 4]) /(cifEST[cifloc,2]#LOG(cifEST[cifloc,2])));
         cifEST=cifEST[, (1:2||4)];
         cifEST = cifEST || lowci || upci;
         create &out.plot from cifEST[colname={"&tname" "CIF" "StdErr" "LowerCI" "UpperCI"}];
         end;
      else do;
         /*Delta Method -- Marubini*/
         lowci[cifloc] = cifEST[cifloc, 2]##exp((-z*cifEST[cifloc, 3]) /(cifEST[cifloc,2]#LOG(cifEST[cifloc,2])));
         upci[cifloc] = cifEST[cifloc, 2]##exp((z*cifEST[cifloc, 3]) /(cifEST[cifloc,2]#LOG(cifEST[cifloc,2])));
         cifEST=cifEST[,(1:3)];
         cifEST = cifEST|| lowci || upci;
         create &out.plot from cifEST[colname={"&tname" "CIF" "StdErr" "LowerCI" "UpperCI"}];
         end;



      append from cifEST;
      close &out.plot;

      create _outgrpLabel from grpLabel[colname={"&gname"}] ;
                append  from grpLabel;
      close _outgrpLabel;

          if (nstrata>1) then do;
                  create outstratLabel from stratLabel[colname={"&stratname"}] ;
                        append  from stratLabel;
                  close outstratLabel;
                  end;


   end;


   /* start Gray's K-sample test, and strata is allowed;
   * merge the data;*/


   if ( ng>1 &(^&notestflag)) then do;
      /* recode the group member id from 1 to ng*/
      grpcode=j(nrow(x), 1, 0);
      do g= 1 to ng;
        grpcode[loc(group= ugg[g])]=g;
      end;



      do k= 1 to nstrata;
         * subset of data for stratum k;
         xsub= x[loc(strat= ust[k]), ];
                 grpcodesub=grpcode[loc(strat= ust[k])];
                 if (ncol(unique(grpcodesub))=ng) then do;

                        Run test_ksample(xsub[,1], xsub[,4], grpcodesub, &rho);
                        score_st= score_st+ score;
                        vcov_st= vcov_st+ vcov;
                        end;
                 else
                        call symput('inconsistentgrpnumflag',left(char(1)));


      end;


      pos= 0;
      do i= 1 to ng1;
         do j= 1 to i;
            pos= pos+ 1;
            vv[i,j]= vcov_st[pos];
            vv[j,i]= vv[i,j];
         end;
      end;

      * test statistic and p-value;

      test_stat= score_st`* ginv(vv)* score_st;
      df= ng1;
      pval= 1- probchi(test_stat, df);
      Test=  test_stat || df || pval;
      create TestResult from Test[colname={ "test_stat" "df" "pval"}];
      append from Test;
      close TestResult;

      end;

quit; /*Exit Proc IML*/


/*create the output data for the cumulative incidence estimates for the CIF plot */

%if &nocifflag %then %do; %put WARNING: No CIF is computed when both STRATA and GROUP are specified; %goto exit2; %end;
data &out.plot;
   merge _outgrpLabel &out.plot;
   if (0) then set &data(keep=&tname);
run;


%if (&nstrata>1) %then %do;
  data &out.plot;
    label &stratname="&stratlabel";
    merge outstratLabel &out.plot;
  run;
%end;

data &out(%if %bquote(&group)= %then drop=&gname;);
   set &out.plot;
   %if (&nstrata>1) %then %do;
      by &stratname &gname;
          %end;
   %else %do;
          by &gname;;
      %end;
   if not last.&gname;


run;

%if( ^&nocifestflag ) %then  %do;
Title "Cumulative Incidence Estimates  with &Confidence% Confidence Limits";
proc print data=&out NOOBS Label;
   label   LowerCI= "Lower &Confidence% Limits"
           UpperCI="Upper  &Confidence% Limits";
run;
%end;

   %if &notestflag ~= 1 AND &group ~= %STR() %then %do;
      /* Get Gray's p-value */
      PROC SQL noprint;
         select pval format=pvalue6.4 into :pval
	     from testresult;
      QUIT;
   %end;

/*plot the CIF curves*/

%if &noplotflag %then %goto exit2;
    %if (%bquote(&title)=) %then
        %if (&plotclflag ) %then
            %let title= Cumulative Incidence Functions with &Confidence% Confidence Intervals;
            %else
                %let title= Cumulative Incidence Functions;;


title '&title';
title2 '&title2';

%if (&ng>1 ) %then %let groupopt=group=&gname;
%else %let groupOpt=;

%if (&nstrata=1) %then %do;
        proc sgplot data=&out.plot;
           label CIF='Cumulative Incidence';
           yaxis values=(0 .2 .4 .6 .8 1.0);
           step y=CIF x=&tname / &groupOpt  name="cif";
		   %if &notestflag ~= 1 and &group ~= %STR() %then %do; 
		      inset "Gray's p=&pval"/position=topright noborder; 
           %end;
           %if &plotclflag %then
                band x=&tname lower=LowerCI upper=UpperCI / group=&gname fill transparency=.5
                                                     modelname='cif' legendlabel='';;

        run;
        %end;
%else %do;

        proc sgpanel data=&out.plot;
           panelby &stratname;
           label CIF='Cumulative Incidence';
           *rowaxis values=(0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0);

           step y=CIF x=&tname / &groupOpt  name="cif";
           %if &plotclflag %then
               band x=&tname lower=LowerCI upper=UpperCI / group=&gname fill transparency=.5
                                                          modelname='cif' legendlabel='';;

        run;
        %end;



%exit2:
   %if (&inconsistentgrpnumflag) %then %do; %put WARNING: The test is not available because the &group number is not the same across the strata; %goto exit; %end;
   %if ( &ng>1 &(^&notestflag) AND &noprinttestflag=0) %then %do;
       Title "Gray's Test for Equality of Cumulative Incidence Functions";
       proc print data=testresult NOOBS Label;
       label test_stat="Chi-Square"
                 df="DF"
                 pval= "Pr > Chi-Square";
           format pval pvalue6.4;
       run;
   %end;

%exit:
Title ;
Title2 ;
proc datasets nolist;
     delete &data.1 _outgrpLabel _c;
run;

options &savenote;


%mend CIF;


