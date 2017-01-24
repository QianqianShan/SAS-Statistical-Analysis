/* Survival Analysis Using SAS: A Practical Guide, Second Edition
By Paul D. Allison */



/* Chapter 3 Estimating and comparing Survival Curves with PROC LIFETEST */

/* The following program is found on page 31 of */
/* "Survival Analysis Using SAS"    */

/*An example of Kaplan Meier method */
DATA myel;
INPUT dur status treat renal;
DATALINES; 
   8 1 1 1
 180 1 2 0
 632 1 2 0
 852 0 1 0
  52 1 1 1
2240 0 2 0
 220 1 1 0
  63 1 1 1
 195 1 2 0
  76 1 2 0
  70 1 2 0
   8 1 1 0
  13 1 2 1
1990 0 2 0
1976 0 1 0
  18 1 2 1
 700 1 2 0
1296 0 1 0
1460 0 1 0
 210 1 2 0
  63 1 1 1
1328 0 1 0
1296 1 2 0
 365 0 1 0
  23 1 2 1
;
PROC LIFETEST DATA=myel method=KM;	*method=KM is default;
	 	TIME dur*status(0);
RUN;
/*status(0) means 0 is the corresponding value of status for censored obs 
  TIME statement: The TIME statement is required. 
                  It is used to indicate the failure time variable,*/



/* The following program is found on page 34 of */
/* "Survival Analysis Using SAS"      */
/* Use output delivery system to produce survival density funtion (plot=s)*/
ODS GRAPHICS ON;
PROC LIFETEST DATA=myel PLOTS=S;
	TIME dur*status(0);
RUN;
ODS GRAPHICS OFF;


/* The following program is found on page 35 of */
/* "Survival Analysis Using SAS"      */
/*options to not show censored data and show confidence intervals*/
ODS GRAPHICS ON;
title 'Plot without showing censored data';
PROC LIFETEST DATA=myel PLOTS=survival(atrisk nocensor cl);
TIME dur*status(0);
run;
title;
ODS GRAPHICS OFF;


/*Two ways to obtain confidence bands and pointwise confident limits:

Hall-wellner(hw) method and equal precision (ep),
ep is more stable at tails */

ODS GRAPHICS ON;
title 'Plot with confidence bands';
PROC LIFETEST DATA=myel PLOTS=s(cl cb=ep);
TIME dur*status(0);
run;
title;
ODS GRAPHICS OFF;
/*cb=all,ep,hw */


/* The following program is found on page 37 of */
/* "Survival Analysis Using SAS"    */

/*write survival estimates and POINTWISE confidence limits into a dataset*/
title;
PROC LIFETEST DATA=myel OUTSURV=a;
	TIME dur*status(0);
PROC PRINT DATA=a;
RUN;
/*outsurv=Names an output data set to contain survival estimates and confi-
dence limits,
confidence bands always wider than pointwise confidence limits*/


/* The following program is found on page 39 of */
/* "Survival Analysis Using SAS"      */
/*Test for difference in survivor functions*/
title 'define strata by strata';
ODS GRAPHICS ON;
PROC LIFETEST DATA=myel PLOTS=S(TEST);
	TIME dur*status(0);
	STRATA treat;
RUN;
ODS GRAPHICS OFF;
/*The STRATA statement identifies the variables that determine the strata levels.
  Strata are formed according to the nonmissing values of these variables*/



title 'define strata by BY';
/*Need first sort the data in ascending order */
proc sort data=myel out=myel1;
by treat;
run;
proc print data=myel1;
run;
ODS GRAPHICS ON;
PROC LIFETEST DATA=myel1 PLOTS=S(TEST);
	TIME dur*status(0);
	BY treat;
RUN;
ODS GRAPHICS OFF;
/*The BY statement is more efficient than the STRATA statement for defining strata 
in large data sets. 
However,if you use the BY statement to define strata, PROC LIFETEST does not pool 
over strata for testing the association of survival time with covariates,
nor does it test for homogeneity across the BY groups.*/



/* The following program is found on page 44 of  */
/* "Survival Analysis Using SAS"    */
title 'define strata by strata and display all tests except for log-rank and wilcoxn';
ODS GRAPHICS ON;
PROC LIFETEST DATA=myel PLOTS=S(TEST);
	TIME dur*status(0);
	STRATA treat/ TESTS=ALL;
RUN;
ODS GRAPHICS OFF;



/* The following program is found on page 46 of  */
/* "Survival Analysis Using SAS"    */
* libname mylib 'C:\Users\Qianqian Shan\Dropbox\';

/*read sas7bdat file*/
data myrecid;
set 'C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis\recid.sas7bdat';
run;
proc print data=myrecid(obs=10);
run;
PROC LIFETEST DATA=myrecid plots=s;
	TIME week*arrest(0);
	STRATA wexp paro / ADJUST=TUKEY;
RUN;
/*adjust: produce p values for six pairwise comparisons fo the four strata dand then 
to report p values that have been adjust for multiple comparison using TUkey's method*/

/* The following program is found on page 47 of  */
/* "Survival Analysis Using SAS"    */

PROC LIFETEST DATA=myrecid plots=s;
	TIME week*arrest(0);
	STRATA age(21 24 28)/adjust=bon;
RUN;
/*For numeric values such as age here, can define groups by intervals rather than
   unique values. each stratum is identified by the midpoint of the interval except 
   for the two ending intervals. */




/* The following program is found on page 49 of  */
/* "Survival Analysis Using SAS"    */
/*When the number of obs is large and event times are precisely measured, 
  many unique values will be generated and KM method will produce long tables 
  which not preferred

Here we used life table which divides the data into several intervals*/
title 'life table';
PROC LIFETEST DATA=myrecid METHOD=LT; *or method=life;
	TIME week*arrest(0);
RUN;



/* The following program is found on page 53 of */
/* "Survival Analysis Using SAS"    */

/*Get survival and hazard estimates by plots=(s,h)*/
title 'plots';
ODS GRAPHICS ON;
PROC LIFETEST DATA=myrecid METHOD=LIFE PLOTS=(S,H);
	TIME week*arrest(0);
RUN;
ODS GRAPHICS OFF;

/*In the above life table, as the data is censored at week 52, but the interval is [50,60),
   introduced some overestimation.*/

/* The following program is found on page 54 of */
/* "Survival Analysis Using SAS"    */

/*Explicitly set the last interval to be [50,53) so week 52 is included. */
title 'fixed interval';
DATA newrecid;
	SET myrecid;
	IF arrest=0 THEN week=53;
PROC LIFETEST DATA=newrecid METHOD=LIFE PLOTS=(S,H)
	INTERVALS=10 20 30 40 50 53;
	TIME week*arrest(0);
RUN;
/*the drop on hazard function is less steep*/


/* The following program is found on pages 56-57 of */
/* "Survival Analysis Using SAS"    */

/*Construct life tables from published data which provides only :
  boundaries of the intervals, number of events,number of censored cases */
title 'Life table from grouped data';
DATA grouped;
INPUT time status number;
DATALINES;
25 1 16
25 0 3
75 1 11
75 0 0
150 1 4
150 0 2
300 1 5
300 0 4
550 1 2
550 0 6
850 1 4
850 0 3
1150 1 1
1150 0 2
1450 1 1
1450 0 3
1700 1 0
1700 0 1
;
ods graphics on;
PROC LIFETEST data=grouped METHOD=LIFE INTERVALS=50 100 200 400 700
   1000 1300 1600 PLOTS=(survival(cl) h);
	TIME time*status(0);
	FREQ number;
RUN;
ods graphics off;
/*FREQ statement is used as a weight variable*/


/* The following program is found on page 60 of */
/* "Survival Analysis Using SAS"    */

title 'Test for effects of covariates';
PROC LIFETEST DATA=myrecid;
	TIME week*arrest(0);
	TEST fin age race wexp mar paro prio;
RUN;
/*The TEST statement specifies a list of numeric covariates (prognostic variables) 
  that you want tested for association with the failure time.
Two sets of rank statistics are computed. These rank statistics and their variances
 are pooled over all strata.
Univariate (marginal) test statistics are displayed for each of the covariates.*/



/* The following programs are found on page 63 of */
/* "Survival Analysis Using SAS"    */

title 'strata argument';
PROC LIFETEST DATA=myel;
	TIME dur*status(0);
	STRATA renal;
	TEST treat;
RUN;

title 'group statement';
PROC LIFETEST DATA=myel;
	TIME dur*status(0);
	STRATA renal / GROUP=treat;
RUN;
/*When there are both STRATA and TEST statements, log-rank and Wilcoxon statistics produced 
   by TEST are first calculated within STRATA and then averaged across strata. 

STRATA and GROUP get same chi-square values, but GROUP can have more than two categories, while 
 TEST will treat variables with more than two values as quantitative measure. */

/* The following program is found on page 66 of */
/* "Survival Analysis Using SAS"    */


/*Log survival and smoothed hazard plots */
title 'Smooth hazard function with ODS';
ODS GRAPHICS ON;
PROC LIFETEST DATA=myrecid PLOTS=H;
	TIME week*arrest(0);
RUN;
ODS GRAPHICS OFF;



/* The following program is found on page 67 of */
/* "Survival Analysis Using SAS"    */
title 'Smooth hazard function with ODS and bandwidth=5 weeks';
proc print data=myrecid(obs=10);
run;
ODS GRAPHICS ON;
PROC LIFETEST DATA=myrecid PLOTS=H(bw=5);
	TIME week*arrest(0);
RUN;
ODS GRAPHICS OFF;



/* The following program is found on page 68 of */
/* "Survival Analysis Using SAS"    */
title 'Smooth hazard function with ODS and CL';
ODS GRAPHICS ON;
PROC LIFETEST DATA=myrecid PLOTS(only)=H(cl);
	TIME week*arrest(0);
RUN;
ODS GRAPHICS OFF;





/*Chapter 4 Estimating parametric regression models with PROC LIFEREG */
/*LIFEREG 1. accomadates left censoring and interval censoring
          2. can test cenrtain hypotheses about the shape of the hazard function
  PHREG  1. allows only right censoring 
         2. only gives non-parameric estimates of the survivor functions, which 
            can be difficult to interpret. 
 distribution= exponenetial gamma llogistic lnormal logistic normal weibull*/

/* The following program is found on page 74 of */
/* "Survival Analysis Using SAS"    */

title 'Regression with lognormal model';
PROC LIFEREG DATA=myrecid;
	MODEL week*arrest(0)=fin age race wexp mar paro prio
		/ DISTRIBUTION=LNORMAL;
RUN;

title 'Regression with Weibull model';
PROC LIFEREG DATA=myrecid;
	MODEL week*arrest(0)=fin age race wexp mar paro prio
		/ DISTRIBUTION=weibull;
RUN;

title 'Regression with gamma model';
PROC LIFEREG DATA=myrecid;
	MODEL week*arrest(0)=fin age race wexp mar paro prio
		/ DISTRIBUTION=gamma;
RUN;
/*week is right censored data */



/* The following program is found on page 87 of */
/* "Survival Analysis Using SAS"    */

/*categorical variables and the CLASS statement */
title 'Regression with educ as categorical variable';
PROC LIFEREG DATA=myrecid;
	CLASS educ;
	MODEL week*arrest(0)=fin age race wexp mar paro
		prio educ / D=WEIBULL;
RUN;



/* The following program is found on page 96 of */
/* "Survival Analysis Using SAS"    */

/*Fit a null model to calculate te statistic when assuming all cavariates' coefficients are zero*/
title 'Null model';
proc LIFEREG data=myrecid; 
MODEL week*arrest(0)= / D=WEIBULL;
probplot;
run;


/* The following program is found on page 98 of */
/* "Survival Analysis Using SAS"    */


/*Test hypothesis by beta_3=beta_4 */

title 'Test if two betas are equal';
data testrecid;
set myrecid;
IF educ = 3 THEN educ = 4;
run;

proc print data=testrecid(obs=10);
run;

proc lifereg data=testrecid;
class educ;
model week*arrest(0)=fin age race wexp mar paro prio educ / D=WEIBULL; 
PROBPLOT;
run;



/* The following programs are found on page 104 of */
/* "Survival Analysis Using SAS"    */


/*Left censoring and right censoring 
MODEL (lower,upper)=list of covariates;
i.e. use two time variables
*/

title 'Use LIFEREG to deal with left and interval censored data';
PROC SORT DATA=myrecid OUT=recid2;
	BY DESCENDING arrest;
DATA recidlft;
	SET recid2;
	IF _N_ LE 30 THEN week = .;
RUN;



DATA recid3;
	SET recidlft;
		/* uncensored cases: */
	IF arrest=1 AND week ne . THEN DO;
		upper=week;
		lower=week;
	END;
		/* left-censored cases: */
	IF arrest=1 AND week = . THEN DO;
		upper=52;
		lower=.;
	END;
		/* right-censored cases: */
	IF arrest=0 THEN DO;
		upper=.;
		lower=52;
	END;
RUN;

PROC LIFEREG DATA=recid3;
	MODEL (lower,upper)=fin age race wexp mar paro prio
		/ D=WEIBULL;
RUN;


/* The following program is found on page 107 of */
/* "Survival Analysis Using SAS"    */

/*As week is not an exact day, can also be considered as an interval-censored data */
title 'Use week data as interval censored';
DATA recidint;
	SET myrecid;
		/* interval-censored cases: */
	IF arrest=1 THEN DO;
		upper=week;
		lower=week-.9999;
	END;
		/* right-censored cases: */
	IF arrest=0 THEN DO;
		upper=.;
		lower=52;
	END;
RUN;
PROC LIFEREG DATA=recidint;
	MODEL (lower, upper) = fin age race wexp mar paro prio
		/ D=WEIBULL;
RUN;



/* The following program is found on pages 108-109 of */
/* "Survival Analysis Using SAS"    */
title 'Generate the predicted median survival time as point estimates for each individual';
PROC LIFEREG DATA=myrecid;
	MODEL week*arrest(0)=fin age race wexp mar paro prio
		/ D=WEIBULL;
	OUTPUT OUT=a P=median STD=s cdf=cdf1;
RUN;

PROC PRINT DATA=a;
	VAR week arrest _prob_ median s cdf1;
RUN;

/*TWO MARCROS USED */
/* The following program is found on pages 294-295 of */
/* "Survival Analysis Using SAS"*/    

%macro lifehaz(outest=,out=,obsno=0,xbeta=lp);

/********************************************************************
Version 2.0 (9-14-01)

This version of LIFEHAZ works for SAS Release 6.12 through 
Release 9.2.  

********************************************************************/

data;
  set &outest;
  call symput('time',_NAME_);
run;
proc means data=&out noprint;
  var &time &xbeta;
  output out=_c_ min(&time)=min max(&time)=max mean(&xbeta)=mean;
run;
data;
  set &outest;
  call symput('model',_dist_);
  s=_scale_;
  d=_shape1_;
  _y_=&obsno;
  set _c_ (keep=min max mean);
  if _y_=0 then m=mean;
  else do;
    set &out (keep=&xbeta) point=_y_;
    m=&xbeta;
  end;
  inc=(max-min)/300;
  g=1/s;
  alph=exp(-m*g);
  _dist_=upcase(_dist_);
if _dist_='LOGNORMAL' or _dist_='LNORMAL'  then do;
  do t=min to max by inc;
  z=(log(t)-m)/s;
  f=exp(-z*z/2)/(t*s*sqrt(2*3.14159));
  Surv=1-probnorm(z);
  h=f/Surv;
  output;
  end;
end;
else if _dist_='GAMMA' then do;
  k=1/(d*d);
  do t=min to max by inc;
  u=(t*exp(-m))**(1/s);
  f=abs(d)*(k*u**d)**k*exp(-k*u**d)/(s*gamma(k)*t);
  Surv=1-probgam(k*u**d,k);
  if d lt 0 then Surv=1-Surv;
  h=f/Surv;
  output;
  end;
end;
else if _dist_='WEIBULL' or _dist_='EXPONENTIAL' or _dist_='EXPONENT'  then do;
  do t=min to max by inc;
  h=g*alph*t**(g-1);
  output;
  end;
end;
else if _dist_='LLOGISTIC' or _dist_='LLOGISTC' then do;
  do t=min to max by inc;
  h=g*alph*t**(g-1)/(1+alph*t**g);
  output;
  end;
end;
else put 'ERROR:DISTRIBUTION NOT FITTED BY LIFEREG';
run;
proc gplot;
  plot h*t / haxis=axis2 vaxis=axis1 vzero;
  symbol1 i=join v=none c=black;
  axis1 label=(f=titalic angle=90 'Hazard');
  axis2 label=(f=titalic justify=c 'time' f=titalic justify=c "&model");
run; quit;
%mend lifehaz;


/* The following program is found on pages 296-297 of */
/* "Survival Analysis Using SAS"  */
  
%macro predict (outest=, out=_last_,xbeta=,time=);

/********************************************************************

Example:  To get 5-year survival probabilities for every individual
in the sample (assuming that actual survival times are measured in 
years);

%predict(outest=a, out=b, xbeta=lp, time=5).

*********************************************************************/

data _pred_;
_p_=1;
set &outest  point=_p_;
set &out;
lp=&xbeta;
t=&time;
gamma=1/_scale_;
alpha=exp(-lp*gamma);
prob=0;
_dist_=upcase(_dist_);
if _dist_='WEIBULL' or _dist_='EXPONENTIAL' or _dist_='EXPONENT' then prob=exp(-alpha*t**gamma);
if _dist_='LOGNORMAL' or _dist_='LNORMAL' then prob=1-probnorm((log(t)-lp)/_scale_);
if _dist_='LLOGISTIC' or _dist_='LLOGISTC' then prob=1/(1+alpha*t**gamma);
if _dist_='GAMMA' then do;
  d=_shape1_;
  k=1/(d*d);
  u=(t*exp(-lp))**gamma;
  prob=1-probgam(k*u**d,k);
  if d lt 0 then prob=1-prob;
  end;
drop lp gamma alpha _dist_ _scale_ intercept
     _shape1_ _model_ _name_ _type_ _status_ _prob_ _lnlike_ d k u;
run;
proc print data=_pred_;
run;
%mend predict;


/* The following program is found on page 110 of */
/* "Survival Analysis Using SAS"    */
title 'Generate the probability of surviving to some specified time and hazard function';
PROC LIFEREG DATA=myrecid OUTEST=a plots=all;
	MODEL week*arrest(0) = fin age race wexp mar paro prio
		/ D=WEIBULL;
	OUTPUT OUT=b XBETA=lp;
RUN;

/*macro to predict the 30-week survival probability*/
%PREDICT(OUTEST=a,OUT=b,XBETA=lp,TIME=30)
/* OUTEST= data set contains parameter estimates and the log likelihood for the model.
   XBETA specifies a variable to contain the computed value of x'b, where x is the 
         covariate vector and b is the vector of parameter estimates.*/
proc print data=a;
run;


/* The following program is found on page 111 of */
/* "Survival Analysis Using SAS"    */

/*Macro to produce hazard function */
%LIFEHAZ(OUTEST=a,OUT=b,XBETA=lp)


/* The following program is found on page 113 of */
/* "Survival Analysis Using SAS"    */

DATA quarter;
	SET myrecid;
	keep week arrest fin age race wexp mar paro prio educ quarter j time event; 
	quarter=CEIL(week/13);
	DO j=1 TO quarter;
		time=13;
		event=0;
	IF j=quarter AND arrest=1 THEN DO;
		event=1;
		time=week-13*(quarter-1);
	END;
		OUTPUT;
	END;
RUN;

proc print data=quarter(obs=20);
run;

/* The following program is found on page 114 of */
/* "Survival Analysis Using SAS"    */

title 'piecewise exponential dist';
PROC LIFEREG DATA=quarter;
	CLASS j;
	MODEL time*event(0)=fin age race wexp mar paro prio j
		/ D=EXPONENTIAL COVB;
RUN;
/*COVB : estimated covariance matrix 
  j: the index variable in the do loop, has values 1,2,3,4. 
  Similar to a CLASS variable for piecewise exponential model */


/* The following program is found on page 118 of */
/* "Survival Analysis Using SAS"    */


/*Bayesian Estimation and Testing 
 Advantage: 1. Inference doesn't reply on large sample approximations and thus could be more 
                accurate for small samples. 
            2. Prior information can be used.
            3. Camparison of non-nested models can be accomplished in a systematic framework. */
ODS HTML;
ODS GRAPHICS ON;
title 'Bayesian Analysis';
PROC LIFEREG DATA=myrecid;
	MODEL week*arrest(0)=fin age race wexp mar paro prio /D=WEIBULL;
	BAYES;  *statement needed for bayesian analysis;
RUN; 
ODS GRAPHICS OFF;
ODS HTML CLOSE;
/*The BAYES statement requests a Bayesian analysis of the regression model 
  by using Gibbs sampling.*/




/* Chapter 5 Estimating Cox Regression Models with PROC PHREG */

/*Advantage of Cox method: 
 1. don't require choice of some particular probability distribution to represent survival times
 2. relatively easy to incorporate time-dependent covariates 
 3. permits a kind of stratified analysis which is effective in controling nuisance parameters. */
/* The following program is found on page 129 of */
/* "Survival Analysis Using SAS"    */


data myrecid;
set 'C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis\recid.sas7bdat';
run;

/*Apply partial likelihood method to the recid data */
title 'Partial likelihood';
PROC PHREG DATA=myrecid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio;
RUN;
/*In the output Model Information part, Ties Handling BRESLOW shows the way to handle tied data.*/



/* The following program is found on page 134 of */
/* "Survival Analysis Using SAS"    */

title 'heart transplant example';
OPTIONS YEARCUTOFF=1900;
DATA stan;
   set 'C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis\stan.sas7bdat';
run;

data stan;
set stan;
INPUT dob mmddyy9. doa mmddyy9. dot mmddyy9. dls mmddyy9.     
         dead surg m1 m2 m3;
   surv1=dls-doa;
   surv2=dls-dot;
   ageaccpt=(doa-dob)/365.25;
   agetrans=(dot-dob)/365.25;
   wait=dot-doa;
   IF dot=. THEN trans=0; ELSE trans=1;
RUN;

proc print data=stan;
run;

/* The following program is found on page 135 of */
/* "Survival Analysis Using SAS"    */

/*surv1 is the number of days the patient survive(or censored) since the patient entered 
  the program */
title 'Test tranplatation raise or lower the hazard of death';
PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=trans surg ageaccpt;
RUN;
/*Note: trans(transplant status) is the time dependent covariate as people who enter 
        the program longer tends to have a higher probability to find the donor. 
 This issue will be dealt later this Chapter. */


/* The following program is found on page 136 of */
/* "Survival Analysis Using SAS"    */

/*Now restrict to only patients who have the tranplantation and check what is affecting 
  their longevity. */
title 'Only patiens who have tranplanted';
PROC PHREG DATA=stan;
   WHERE trans=1;
   MODEL surv2*dead(0)=surg m1 m2 m3 agetrans wait dot;
RUN;




/* The following program is found on page 145 of */
/* "Survival Analysis Using SAS"    */

title 'ties=breslow(default)';
PROC PHREG DATA=myrecid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio 
         / TIES=BRESLOW;
RUN;

title 'ties=DISCRETE';
PROC PHREG DATA=myrecid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio 
         / TIES=discrete;
RUN;

title 'ties=efron';
PROC PHREG DATA=myrecid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio 
         / TIES=efron;
RUN;

title 'ties=exact';
PROC PHREG DATA=myrecid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio 
         / TIES=EXACT;
RUN;
/*The EXACT method can take a considerable amount of computer resources. If ties are not
extensive, the EFRON and BRESLOW methods provide satisfactory approximations to the
EXACT method for the continuous time-scale model.

In general, Efron’s approximation gives results that are much closer to the EXACT method 
results than Breslow’s approximation does. If the time scale is genuinely discrete, 
you should use the DISCRETE method. 

The DISCRETE method is also required in the analysis of case-control studies when there is
more than one case in a matched set. 

If there are no ties, all four methods result in the same likelihood and yield
identical estimates. 

The default, TIES=BRESLOW, is the most efficient method when there are no ties

*/




/* The following program is found on page 155 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg ageaccpt / TIES=EFRON;
   IF wait>=surv1 OR wait=. THEN plant=0; ELSE plant=1;
   * i.e. if the patient dies before get tranplanted or not transplanted, plant=0
     The IF statement must follow the MODEL statement;
RUN;


/* The following program is found on page 156 of */
/* "Survival Analysis Using SAS"    */

/*Redo the analysis using the counting process. Need to construct multiple records per person,
  one record for each period of time such that the covariates are constant. 

  For example, the first period is the time between acceptance to date of transplant, the 
  second period is the date of tranplant to either death or cencored date. 

  For patients who didn't have transplantation (trans=0),death2 indicator is equal to the 
   original indicator. stop variable equals to surv1(date last seen - date of acceptance),
   start variable is zero .

  For patients who had transplantation, in the first period, start=0,plant=0,stop is the
   waiting time,dead2=0; in the second period, start is the end of last period, and stop 
   is surv1(the time the patient is either dead or censored).

  If stop=0(patient died when entering the program), SAS will delete the data which have the
   same start and stop time, so we need to change the stop to a small nonzero variable. 
 */
DATA stanlong;
SET stan;
plant=0;
start=0;
IF trans=0 THEN DO;
  dead2=dead;
  stop=surv1;
  IF stop=0 THEN stop=.1;
  OUTPUT; 
END;
ELSE DO;
  stop=wait;
  IF stop=0 THEN stop=.1;
  dead2=0;
  OUTPUT;
  plant=1;
  start=wait;
  IF stop=.1 THEN start=.1;
  stop=surv1;
  dead2=dead;
  OUTPUT;
END;
RUN;
/*As a result, stanlong will have one record for patiens who did not tranplant, two for trans*/


/* The following program is found on page 158 of */
/* "Survival Analysis Using SAS"    */

title 'Regression using counting process';
PROC PHREG DATA=stanlong;
   MODEL (start,stop)*dead2(0)=plant surg ageaccpt /TIES=EFRON;
RUN;
/*Note results are the same as :

PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg ageaccpt / TIES=EFRON;
   IF wait>=surv1 OR wait=. THEN plant=0; ELSE plant=1;
RUN;

Note: this is is different with the interval censored data in LIFEREG process. 
      In LIFEREG, the event is known to occur sometime within the specified interval,
      While in counting process syntax, (start,stop) represents the interval of time 
        during which the individual was continuously at risk of event and events can 
        ONLY occur at the end of the interval. 
*/


/* The following program is found on page 160 of */
/* "Survival Analysis Using SAS"    */


/*try to use an alternative age with different time origin*/
title 'Regression with different current age';
PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg age / TIES=EFRON;
   IF wait>surv1 OR wait=. THEN plant=0; ELSE plant=1;
   age=ageaccpt+surv1;
RUN;
/*The results have no difference with the original age we used
  due to the reason that we are assuming log(h(t)) is a linear function of age.*/


/* The following program is found on page 161 of */
/* "Survival Analysis Using SAS"    */

title 'Regression with log(age)';
PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg logage / TIES=EFRON;
   IF wait>surv1 OR wait=. THEN plant=0; ELSE plant=1;
   logage=LOG(ageaccpt+surv1);
RUN;
title;
/*logage=log(age)*/


/* The following programs are found on page 162 of */
/* "Survival Analysis Using SAS"    */
* filename rawinput "C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis\recid.sas7bdat";
data recid;
set "C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis\recid.sas7bdat";
keep week arrest fin age race wexp mar paro prio educ emp1-emp52;
run;
/*emp1-emp51 are is coded as 1 is the the person is fulled employed at the 
   corresponding week */
proc print data=recid(obs=10);
run;

title;
PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   DO i=1 TO 52;
      IF week=i THEN employed=emp[i];
   END;
RUN;
/* ARRAY statement treats the 52 distinct dummy variables as a single subscripted array.
   Another way to create employed indicator column is shown below. */





/* The following program is found on page 163 of */

/*Another way to create employed variable*/
PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employed=emp[week];
RUN;



/* The following programs are found on page 164 of */
/* "Survival Analysis Using SAS"    */

/*Use employment of the previous week for lgging covariate values */
PROC PHREG;
   WHERE week>1;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employed=emp[week-1];
RUN;
/*A little better model*/

/*A model with 2 lagged weeks */
PROC PHREG;
   WHERE week>2;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employ1 
         employ2 / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employ1=emp[week-1];
   employ2=emp[week-2];
RUN;
/*Neither employ1 or employ2 are significant as they are highly correlated. */



/* The following programs are found on page 165 of */
/* "Survival Analysis Using SAS"    */

/*Check the possibility that the hazard of arrest may depend on the cumulative 
 employment experience */
DATA recidcum;
   SET recid;
   ARRAY emp(*) emp1-emp52;
   ARRAY cum(*) cum1-cum52;
   cum1=emp1;
   DO i=2 TO 52;
      cum(i)=(cum(i-1)*(i-1) + emp(i))/i;
	  *cum is a proportion;
   END;

   proc contents data=recid;
   run;

   proc contents data=recidcum;
   run;

   proc print data=recidcum;
   run;

PROC PHREG DATA=recidcum;
   WHERE week>1;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employ 
         / TIES=EFRON;
   ARRAY cumemp(*) cum1-cum52;
   EMPLOY=cumemp[week-1];
RUN;

/* Counting process method */
DATA recidcount;
  ARRAY emp(*) emp1-emp52;
  SET recid;
  arrest2=0;
  DO stop=1 TO week;
    start=stop-1;
    IF stop=week THEN arrest2=arrest;
    employed=emp(stop);
    *IF week>=1 THEN emplag=emp(stop-1);
     *ELSE emplag=.;
    OUTPUT;
  END;
RUN;

proc contents data=recidcount;
run;

proc print data=recidcount(obs=100);
run;


/* The following program is found on page 166 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=recidcount;
MODEL (start,stop)*arrest2(0)=fin age race wexp mar paro prio employed /TIES=EFRON;
RUN;
/*This produces the same results as 

PROC PHREG;
   WHERE week>1;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employed=emp[week-1];
RUN;

*/




/* The following program is found on page 167 of */
/* "Survival Analysis Using SAS"    */


/*Ad-hoc estimates of time dependent covariates 

 We may have covariates measured at regular intervals but the intervals don't corresponds 
  to the units in which event times are measured. For example, we know the exact day of 
  failure but don't the one covariate at that time as it's measured only monthly. */
DATA blood;
   set "C:\Usersblood.dat"; *no such data available;
   keep deathday status alb1-alb12;
run; 
PROC PHREG;
   MODEL deathday*status(0)=albumin;
   ARRAY alb(*) alb1-alb12;
   deathmon=CEIL(deathday/30.4);
   albumin=alb[deathmon];
RUN;



/* The following program is found on pages 167-168 of */
/* "Survival Analysis Using SAS"    */

DATA bloodcount;
SET blood;
  ARRAY alb(*) alb1-alb12;
status2=0;
  deathmon=CEIL(deathday/30.4);
DO j=1 TO deathmon;
  start=(j-1)*30.4;
  stop=start+30.4;
  albumin=alb(j);
  IF j=deathmon THEN DO;
    status2=status;
    stop=deathday-start;
  END;
  OUTPUT;
END;
PROC PHREG DATA=bloodcount;
   MODEL (start,stop)*status2(0)=albumin;
RUN;



/* The following programs are found on page 168 of */
/* "Survival Analysis Using SAS"    */
PROC PHREG DATA=blood;
   MODEL deathday*status(0)=albumin;
   ARRAY alb(*) alb1-alb12;
   deathmon=deathday/30.4;
   j=CEIL(deathmon);
   IF j=1 THEN albumin=alb(1);
   ELSE albumin=alb[j]+(alb[j]-alb[j-1])*(deathmon-j+1);
RUN;

DATA bloodcount;
SET blood;
  ARRAY alb(*) alb1-alb12;
status2=0;
DO i=1 TO deathday;
  start=i-1
  stop=i;
  j=CEIL(i/30.4);
  albumin=alb(j)+(alb(j)-alb(j-1))*(i/30.4-j+1);
  IF i=deathday THEN status2=status;
  OUTPUT;
END;
PROC PHREG DATA=bloodcount;
   MODEL (start,stop)*status2(0)=albumin;
RUN;



/* The following program is found on pages 170-171 of */
/* "Survival Analysis Using SAS"*/


/*Time dependent covariates that change at irregular intervals */
/*Data explanation: 
1. surv is either the time of death or time of censoring; 
2. time1-tim10: the time to visit the clinic 
3. pt1-pt10: coagulation time(PT) level measure each visit. 
*/
libname mylib 'C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis';
PROC PHREG DATA=mylib.alco;
   MODEL surv*dead(0)=pt;
   time1=0;
   ARRAY time(*) time1-time10;
   ARRAY p(*) pt1-pt10;
   DO j=1 TO 10;
   IF surv > time[j] AND time[j] NE . THEN pt=p[j];
   END;
RUN;
/*Not significant pt effects. */
data alco1;
set mylib.alco;
time1=0;
ARRAY time(*) time1-time10;
ARRAY p(*) pt1-pt10;
DO j=1 TO 10;
   IF surv > time[j] AND time[j] NE . THEN pt=p[j];
   END;
run;
proc print data=alco1;
run;
/*Note : 
 In the DO END loop, pt will be reassigned by p[i] until the conditions are not satisfied. 

It's NOT generating multiple data. */


/*Now try to check the relationship between survival time and the relative CHANGE of pt level*/
PROC PHREG DATA=mylib.alco;
   MODEL surv*dead(0)=pt;
   time1=0;
   ARRAY time(*) time1-time10;
   ARRAY p(*) pt1-pt10;
   DO j=1 TO 10;
   IF surv > time[j] AND time[j] NE . THEN pt=p[j]-pt1;
   END;
RUN;
/* Significant */


/* The following program is found on page 171 of */
/* "Survival Analysis Using SAS"    */


/*Realize the above process by counting process. */
DATA alcocount;
 SET mylib.alco;
 time1=0;
 time11=.;
 ARRAY t(*) time1-time11;
 ARRAY p(*) pt1-pt10;
 dead2=0;
 DO j=1 TO 10 WHILE (t(j) NE .);
   start=t(j);
   pt=p(j)-pt1;
   stop=t(j+1);
   IF t(j+1)=. THEN DO;
      stop=surv;
	dead2=dead;
   END;
   OUTPUT;
END;
keep surv start stop dead2 pt;
run;

proc print data=alcocount;
run;


PROC PHREG DATA=alcocount;
  MODEL (start,stop)*dead2(0)=pt;
RUN;
/*Results are the same */



/* The following program is found on page 173 of */
/* "Survival Analysis Using SAS"    */


/*Cox Models with Non-proportional hazards. 

Test the proportional hazard assumption by ASSESS */
ODS GRAPHICS ON;
PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio
         / TIES=EFRON;
   ASSESS PH / RESAMPLE;
RUN;
ODS GRAPHICS OFF;
/* ph is proportional hazard 
 RESAMPLE requests that the Kolmogorov-type supremum test be computed on 1,000 
  simulated patterns or on n simulated patterns if n is specified

The plots generated by ASSESS: the solid line is the observed empirical score process,
 the dashed lines are empirical score process based on 20 random simulations that are 
 based on the assumption of proportional hazards. */





/* The following program is found on page 176 of */
/* "Survival Analysis Using SAS"    */


/*Schoenfeld residual is used to detect possible departures from the proportional 
   hazards. Test if this residual is correlated with time is equivalent as testing
   for ph assumption. 

First save the residual into a file called 'b' */
PROC PHREG DATA=recid;
  MODEL week*arrest(0)=fin age prio race wexp mar paro prio /
    TIES=EFRON;
OUTPUT OUT=b RESSCH=schfin schage schprio schrace schwexp schmar 
  schparo schprio;
RUN;
/*RESSCH specifies the Schoenfeld residuals. These residuals are useful in assessing
   the proportional hazards assumption.
the variables after RESSCH are arbitrarily names corresponding to the eight covariates.*/

proc print data=b(obs=20);
run;


/*Modify the data to include two familiar functions of time: log and the square */
DATA c;
  SET b;
  lweek=log(week);
  week2=week**2;
PROC CORR;
VAR week lweek week2 schfin schage schprio schrace schwexp schmar schparo;
RUN;
/*Results are consistent, however, note that the sample size is only 114 
  (S residual not defined for censored data). */


/* The following program is found on page 178 of */
/* "Survival Analysis Using SAS"    */

/*Interactions with time as time-dependent covariates. 

The following code shows how to estimate Cox model with the effect of age change with time*/

PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio ageweek
         / TIES=EFRON;
   ageweek=age*week;
RUN;
/*
The model: log(h(t))=alpha(t)+(beta_1+beta_2*t)x
ASSESS is not available for time dependent covariates */

/*Test the effect of age is significantly different from 0 at any given week=30 */ 
data recid;
set mylib.recid;
run;

PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio ageweek
         / TIES=EFRON;
   ageweek=age*week;
   test age+30*ageweek=0;
RUN;

proc print data=recid(obs=3);
run;

/*Summary: approches for non-proportionality 
 1. add a time dependent covariate representing the interaction of the original covariate
     and time 
 2. stratification .

*/

/* The following program is found on page 180 of */
/* "Survival Analysis Using SAS"    */


DATA myel;
INPUT dur status treat renal;
DATALINES; 
   8 1 1 1
 180 1 2 0
 632 1 2 0
 852 0 1 0
  52 1 1 1
2240 0 2 0
 220 1 1 0
  63 1 1 1
 195 1 2 0
  76 1 2 0
  70 1 2 0
   8 1 1 0
  13 1 2 1
1990 0 2 0
1976 0 1 0
  18 1 2 1
 700 1 2 0
1296 0 1 0
1460 0 1 0
 210 1 2 0
  63 1 1 1
1328 0 1 0
1296 1 2 0
 365 0 1 0
  23 1 2 1
;
run;

data mylib.myel;
set myel;
run;

/* Stratification: ia another approach to non-proportionality. 
It's most useful when the covariate that interacts with time is both categorical 
 and not of direct interest. */

/*Partial likelihood estimation of the model with STATRA*/
PROC PHREG DATA=myel;
   MODEL dur*status(0)=treat / TIES=EXACT;
   STRATA renal;
RUN;



/* The following programs are found on page 184 of */
/* "Survival Analysis Using SAS"    */



/*Left truncation and Late entry into the risk set */

/*agels: age last seen is the age in years last seen (dead or left)*/
DATA stan2;
   set mylib.stan;
   agels=(dls-dob)/365.25;
RUN;


PROC PHREG DATA=stan2;
   MODEL (ageaccpt,agels)*dead(0)=surg ageaccpt / TIES=EFRON;
RUN;
/*Counting process syntax.

  By this design, patients will have no observable risk of death before the acceptance
   (left truncation) as the start time is ageaccpt. 

Note: ageaccpt is both the start time and a covariate */



/* The following program is found on page 185 of */
/* "Survival Analysis Using SAS"    */

/*Add plant as a time dependent indicator of tranplant status */
PROC PHREG DATA=stan2;
   MODEL (ageaccpt,agels)*dead(0)=plant surg ageaccpt / TIES=EFRON;
   IF agetrans>=agels OR agetrans=. THEN plant=0;
   ELSE plant=1;
RUN;


/* The following program is found on page 186 of */
/* "Survival Analysis Using SAS"    */


/*Estimating survivor function */
ODS GRAPHICS ON;
PROC PHREG DATA=mylib.recid PLOTS=S;
   MODEL week*arrest(0)=fin age prio 
         / TIES=EFRON;
   BASELINE OUT=a SURVIVAL=s LOWER=lcl UPPER=ucl;
RUN;
ODS GRAPHICS OFF;
/* The BASELINE statement creates a new SAS data set that contains the baseline function 
   estimates at the event times of each stratum for every set of covariates (x) given in 
   the COVARIATES= data set. If the COVARIATES= data set is not specified, a reference 
   set of covariates consisting of the reference levels for the CLASS variables and 
   the average values for the continuous variables is used. No BASELINE data set is 
   created if the model contains a time-dependent variable defined by means of 
   programming statement.

Survival specifies the survivor function estimate. 
Upper Lower specify the upper and lower confidence limit for survivor function.

PLOTS=s will plot the adjusted survivor function evaluated at the means of covariates. 

Note: the survivor function obtained here is different with the one obtained in LIFETEST
 by Kaplan-Meier, as the former allows for heterogeneity in the hazard/survivor functions while 
 K-M not. */



/* The following program is found on page 189 of */
/* "Survival Analysis Using SAS"    */

ODS GRAPHICS ON;
PROC PHREG DATA=mylib.recid PLOTS(OVERLAY=ROW)=S;
   MODEL week*arrest(0)=fin age prio 
         / TIES=EFRON;
   STRATA fin;
RUN;
ODS GRAPHICS OFF;
/*With BASELINE and STRATA, one can produce graphs of survivor functions in different stratas
   in one plot. */

/* The following program is found on pages 190-191 of */
/* "Survival Analysis Using SAS"    */

/* Use BASELINE to obtain predictions about survial time for a particular set of covariates */
DATA covals;
   INPUT fin age prio;
   DATALINES;
0 40 3
;
PROC PHREG DATA=mylib.recid;
   MODEL week*arrest(0)=fin age prio / TIES=EFRON;
   BASELINE OUT=a COVARIATES=covals SURVIVAL=s LOWER=lcl  UPPER=ucl;
PROC PRINT DATA=a;
RUN;
/*COVARIATES=SAS-data-set
 names the SAS data set that contains the sets of explanatory variable values for which 
 the quantities of interest are estimated. */

/* The following program is found on page 193 of */
/* "Survival Analysis Using SAS"    */


/*Testing linear hypothesis with CONTRAST or TEST */
PROC PHREG DATA=mylib.recid;
   CLASS educ;
   MODEL week*arrest(0)=fin age prio educ/ TIES=EFRON;
   TEST educ3=educ5;
   CONTRAST 'ed3 vs. ed5' educ 0 1 0 -1 ;
RUN;
/*TEST and CONTRAST can do the same thing. 

  CLASS statement will create a set of dummy variables to represent four out of five 
   categories. */

/*By default, the reference level for class is the one with highest value(educ=6),
   can specify it by: */
CLASS educ(REF='2');



/* The following program is found on page 195 of */
/* "Survival Analysis Using SAS"    */


/*Customized hazard ratios */
PROC PHREG DATA=mylib.recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio;
   HAZARDRATIO age / UNITS=5 10;
RUN;
/* Check the hazard ratio for a 5,10 year increase in age. */

/*The HAZARDRATIO statement enables you to request hazard ratios for any variable in 
  the model at customized settings. For example, if the model contains the interaction 
  of a CLASS variable A and a continuous variable X, the following specification 
  displays a table of hazard ratios comparing the hazards of each pair of levels of A
  at X=3: hazardratio A / at (X=3) diff=ALL;*/


/* The following programs are found on page 196 of */
/* "Survival Analysis Using SAS"    */


/*PHREG allows for interactions to be directly specified in the MODEL statment, 
   however, the hazardratio won't be reported unless HAZARDRATIO is specified. */
PROC PHREG DATA=mylib.recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio fin*age;
   HAZARDRATIO fin / at (age=20 25 30 35 40) CL=PL;
RUN;
/*CL=PL requests confidence limits to be computed by profile likelihood, which is 
  more accurate for small samples than conventional Wald intervals. */




/* The following program is found on page 197 of */
/* "Survival Analysis Using SAS"    */

/*Bayes  estimation and testing */
PROC PHREG DATA=mylib.recid;
  MODEL week*arrest(0)=fin age race wexp mar paro prio;
  BAYES;
RUN;


/* Estimate a piecewise exponential model (a special case of Cox model)
   log(h(t))=alpha_j+beta*x_i

 Piecewise exponential model assumas the hazard is constant within each interval but can vary
  across intervals. The time intervals are divided by cut points alpha1, ...
*/
PROC PHREG DATA=mylib.recid;
  MODEL week*arrest(0)=fin age race wexp mar paro prio;
  BAYES piecewise;
RUN;

* BAYES piecewise NBI=0 NMC=1;


/*
Chapter 6  Competing risks 

The characteristic of competing risk is that : the occurence of one type of event removes
 the individual from the risk of all the other event types. 
*/


/* The following program is found on page 209 of */
/* "Survival Analysis Using SAS"    */

/*Estimates and Tests without Covariates */
DATA const;
  SET mylib.leaders;
  event=(lost=1);
  type=1;
DATA nat;
  SET mylib.leaders;
  event=(lost=2);
  type=2;
DATA noncon;
  SET mylib.leaders;
  event=(lost=3);
  type=3;
/*The above data defines three different types of events, constitutional means, natural,
   non-constitutional. */
DATA combine;
  SET const nat noncon;
PROC LIFETEST DATA=COMBINE PLOTS=LLS;
  TIME years*event(0);
  STRATA type;
RUN;
/*PLOTS=LLS requests the log-log survivor function */


/* The following program is found on page 210 of */
/* "Survival Analysis Using SAS"    */


/*Examine smoothed hazard plots using the kernel smoothing option */
ODS GRAPHICS ON;
PROC LIFETEST DATA=combine PLOTS=H(BW=10);
  TIME years*event(0);
  STRATA type;
RUN;
ODS GRAPHICS OFF;


/* The following program is found on page 212 of */
/* "Survival Analysis Using SAS"    */

/*Test for proportionality of hazard functions for different types, 
  equivalent to testing the betas are equal for log(h(t))= .. model */
PROC LOGISTIC DATA=mylib.leaders;
   WHERE lost NE 0; *lost=0 means still in power;
   MODEL lost=years / LINK=GLOGIT;
RUN;
/*LINK=GLOGIT specifies an unordered multinomial logit model instead of the default 
  cumulative logit model. 

 Check the Global Null hypothesis testing and type 3 effects analysis for test of beta=0
  or beta effects. */


/* The following program is found on page 213 of */
/* "Survival Analysis Using SAS"    */


/*Covariate effects via Cox model */
data leaders;
set mylib.leaders;
run;

/*Treat all event types as the same */
PROC PHREG DATA=leaders;
   MODEL years*lost(0)=manner start military age conflict 	
         loginc growth pop land literacy;
   STRATA region;

/* Focus on type 3 and treat treat type 1, 2 as censoring */
PROC PHREG DATA=leaders;
   MODEL years*lost(0,1,2)=manner start military age conflict 	
         loginc growth pop land literacy;
   STRATA region;

/*Focus on type 2 and treat type 1 3 as censoring */
PROC PHREG DATA=leaders;
   MODEL years*lost(0,1,3)=manner start military age conflict 
         loginc growth pop land literacy;
   STRATA region;

/*Focus on type 1 */
PROC PHREG DATA=leaders;
   MODEL years*lost(0,2,3)=manner start military age conflict 	
         loginc growth pop land literacy;
   STRATA REGION;
RUN;
/*From the results, we can see the coefficients differ greatly across different evetn types. 

  Need to do a test for null: beta_j=beta for all types to see if the difference in 
   coefficients is a result of random variation. 

As all the three models use TIES=BRESLOW to handle ties by default, we can add the three 
 -2log() together and minus the -2log() for all types combined to obtain the test statistic
 and the df.

However, if there are ties or we used different methods to handle ties(for example, EFRON),
 the above method won't work, we need to conduct the following test for combined data and 
  -2log(). THe reason to do so is that in this way, the underlying hazard function could be 
  different for the combined data, thus being consistent with the three separate tests.  
   */

/* The following program is found on page 218 of */
/* "Survival Analysis Using SAS"    */


PROC PHREG DATA=combine;
   MODEL years*event(0)=manner start military age conflict 	
         loginc growth pop land literacy / TIES=EFRON;
   STRATA region type;
RUN;


/* The following programs are found on page 221 of */
/* "Survival Analysis Using SAS"    */

/*PHREG only concerns about the rank order of the time variable, so 0's are ok
  LIFEREG excludes any obs wiht times of 0 or less, so we should treat it as left 
  censored data when doing analysis by LIFEREG */
DATA leaders2;
   SET leaders;
   lower=years; 
   upper=years;
   IF years=0 THEN DO;
      lower=.; 
      upper=1; 
   END;
   IF lost IN (0,1,2) THEN upper=.;
RUN;
/*For left censored data, lower=., upper=1, for right censored data upper=. */

PROC LIFEREG DATA=leaders2; 
  CLASS region;
   MODEL (lower,upper)= manner start military age conflict 	
         loginc literacy region / D=EXPONENTIAL;
RUN;    


/* The following program is found on page 225 of */
/* "Survival Analysis Using SAS"    */

DATA leaders3;
   SET leaders;
   lyears=LOG(years+.5); *avoid log(0);
PROC LOGISTIC DATA=leaders3;
   WHERE lost=1 OR lost=3;
   CLASS region / PARAM=GLM;
   MODEL lost=lyears manner age start military conflict loginc
         literacy region;
RUN;
/*WHERE statement onlu selects lost=1 and lost=3 types for analysis. 

  Check "Testing for Global null hypothesis " for the test of all beta's are zero */


/* The following program is found on page 228 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=leaders;	
   MODEL years*lost(0,2)=manner age start military conflict   
      loginc literacy / TIES=EXACT;
   STRATA region;
RUN;


/* The following program is found on page 230 of */
/* "Survival Analysis Using SAS"    */

/*A macro to call cumulative incidence function */
%CUMINCID(DATA=leaders,
           TIME=years,
           STATUS=lost,
           EVENT=3,
           COMPETE=1 2,
           CENSORED=0)






/*
Chapter 7 Analysis of Tied or Discrete Data with LOGISTIC 
*/


/* The following program is found on page 237 of */
/* "Survival Analysis Using SAS"    */
libname mylib 'C:\Users\Qianqian Shan\Dropbox\SAS-Survival-Analysis';
DATA jobyrs;
   SET mylib.jobdur;
   DO year=1 TO dur;
      IF year=dur AND event=1 THEN quit=1; 
      ELSE quit=0;
      OUTPUT;
   END;
RUN;

proc print data=jobyrs;
run;


/* The following program is found on page 238 of */
/* "Survival Analysis Using SAS"    */
 
PROC LOGISTIC DATA=jobyrs;
  CLASS year / PARAM=GLM;
  MODEL quit(DESC)=ed prestige salary year;
RUN;
/*CLASS statement tells that year is a CATEGORICAL variable rather than quantitative. 
  PARAM=glm overrides the default effect coding 
  DESC specifies that the model will predict the probability of a 1 rather than 0 
*/


/* The following program is found on page 244 of */
/* "Survival Analysis Using SAS"    */
 

/*Data with time dependent covariates */
DATA rankyrs;
   set mylib.rank;
   ARRAY arts(*) art1-art10;
   ARRAY cits(*) cit1-cit10;
   IF jobtime=. THEN jobtime=11;
   DO year=1 TO dur;
      IF year=dur THEN promo=event; 
         ELSE promo=0;
      IF year GE jobtime THEN prestige=prest2; 
         ELSE prestige=prest1;
      articles=arts(year);
      citation=cits(year);
      OUTPUT;
   END;
RUN;


proc print data=rankyrs;
run;

proc contents data=rankyrs;
run;
/* The following program is found on page 245 of */
/* "Survival Analysis Using SAS"    */

/*Specify a model with quadratic effects of year */
PROC LOGISTIC DATA=rankyrs;
   MODEL promo(DESC)=undgrad phdmed phdprest articles citation 
         prestige year year*year;
RUN;



/* The following programs are found on page 251 of */
/* "Survival Analysis Using SAS"    */

/*Specify the outcome to be event , 1 for quit, 2 for fired , 0 for neither */
DATA jobyrs2;
   SET mylib.jobdur;
   DO year=1 TO dur;
     IF year=dur THEN outcome=event;
     ELSE outcome=0;
     OUTPUT;
   END;
RUN;

proc contents data=jobyrs2;
run;

proc print data=jobyrs2;
run;

/* Fit a multinomial logit model */
PROC LOGISTIC DATA=jobyrs2;
  MODEL outcome(REF='0')=ed prestige salary year 
    / LINK=GLOGIT;
RUN;


/* The following programs are found on page 253 of */
/* "Survival Analysis Using SAS"    */

/* Do analysis for each different event type separately 

   One can simply eliminate the sample which doesn't belong to the certain type, and then 
     do binomial logit analysis for each type. */
PROC LOGISTIC DATA=jobyrs2;
   WHERE outcome NE 2;
   MODEL outcome(DESC)=ed prestige salary year;
RUN;

PROC LOGISTIC DATA=jobyrs2;
   WHERE outcome NE 1;
   MODEL outcome(DESC)=ed prestige salary year;
RUN;
/*Note: this procedure is justified as a form of conditional maximum likelihood. 

The results are consistent and asympototically normal. 

In practice, both the coefficients and the sd are trivially different from those produced 
  by the simultaneous estimation precedure. 

Advantages for separating the estimation process: 
1. focus on the event types interested 
2. specify quite different models for different event types. */

/*Summmary: 

Maxmimum likelihood methods are used in chapter 7, which is an alternative when there are 
 many ties and many time-dependent covariates. 

*/



/*Chapter 8 Heterogeneity, Repeated Events, and other topics */

/* The following program is found on page 262 of */
/* "Survival Analysis Using SAS"    */


/*Two approaches to deal with repeated data: 
1. do separate analysis for each successive event 
2. treat each interval as a distinct observation, pooling all the intervals together and 
   estimating a single model. 

Problems of approach 2: When two or more obs come from the same unit, they tend to be more 
                        alike than randomly chosen obs. This may lead to standard errors 
                        estimates are biased downward and test statistics biased upward.
Solutions:
2.1 Ignore possible bias and concentrate on better standard error estimates and test statistics
2.2 formulate a model that incorporates unobserved heterogeneity and estimate the model by 
    maximum likelihood or conditional likelihood. 

*/

/*Approach 1 */
data jobmult;
set mylib.jobmult;
run; 

PROC SORT DATA=jobmult;
   BY seq;
proc print data=jobmult(obs=50);
PROC PHREG DATA=jobmult;
   MODEL duration*event(0)=prestige logsal ed;
   BY seq;
RUN;
/*First sort the data by SEQ , so can run separate models for each job in the sequence by
  adding BY statement 
*/


/* The following program is found on page 265 of */
/* "Survival Analysis Using SAS"    */

/*Approach 2 
detect the dependency */
PROC SORT DATA=jobmult;
  BY id seq;
proc print data=jobmult(obs=60);
run;

/*detect the dependency by adding lag1 dur, 
  estimate only the second interval with the first interval as covariate*/
DATA joblag;
  SET jobmult;
  durlag=LAG1(duration);
PROC PHREG DATA=joblag;
  WHERE seq = 2;
  MODEL duration*event(0)=prestige logsal ed durlag;
RUN;
/*DURLAG effect is significant */


/* The following program is found on page 267 of */
/* "Survival Analysis Using SAS"    */

/*Solution 2.1 Robust standard errors 

 COVSANDWICH can correct the dependency when there are repeated obs. 

COVS < (AGGREGATE) >
requests the robust sandwich estimate of Lin and Wei (1989) for the covariance matrix. When
this option is specified, this robust sandwich estimate is used in the Wald tests for testing the
global null hypothesis, null hypotheses of individual parameters, and the hypotheses in the
CONTRAST and TEST statements. In addition, a modified score test is computed in the
testing of the global null hypothesis, and the parameter estimates table has an additional
StdErrRatio column, which contains the ratios of the robust estimate of the standard error
relative to the corresponding model-based estimate. Optionally, you can specify the keyword
AGGREGATE enclosed in parentheses after the COVSANDWICH (or COVS) option, which
requests a summing up of the score residuals for each distinct ID pattern in the computation
of the robust sandwich covariance estimate. This AGGREGATE option has no effects if the
ID statement is not specified.

Advantage: no assumptions about the structure of the dependency. 
However, there is also no correction for bias in the coef that arise from unobserved 
 heterogeneity. 
*/
PROC PHREG DATA=jobmult COVSANDWICH(AGGREGATE);
   MODEL duration*event(0)=prestige logsal ed seq;
   ID id; *needed;
RUN;


/* The following programs are found on page 269 of */
/* "Survival Analysis Using SAS"    */

/*Solution 2.1 coefficients and test statistics */
/*Idea: get coeficients for each successive event and calculate a weighted average of thos 
        coefs using robust covariance matirx for optimal weights. */

/*first create interactions between each covariate and a set of dummy variables which indicate
  whether it's the first event, second...
  then first the model with all these new covariates */
DATA job;
  SET mylib.jobmult;
  ARRAY p (*) p1-p10;
  ARRAY s (*) s1-s10;
  ARRAY e (*) e1-e10;
  DO i=1 TO 10;
    p(i)=prestige*(seq=i);
    s(i)=logsal*(seq=i);
    e(i)=ed*(seq=i);
  END;
RUN;

PROC PHREG DATA=job COVS(AGGREGATE);
   MODEL duration*event(0)=p1-p10 s1-s10 e1-e10 seq;	
   ID id;
   Prestige: TEST p1,p2,p3,p4,p5,p6,p7,p8,p9,p10 / AVERAGE;
   Salary: TEST s1,s2,s3,s4,s5,s6,s7,s8,s9,s10 / AVERAGE;
   Education: TEST e1,e2,e3,e4,e5,e6,e7,e8,e9,e10 / AVERAGE;
RUN;
/*The last three TEST statements return weighted estimates which we need 
COVS() option adjusts standard errors */

/* The following program is found on page 271 of */
/* "Survival Analysis Using SAS"    */


/*Solution 2.2 formulate model to estimate the unobserved heterogeneity */
/*
  Check page 90-93 for details of the math. Here we assume the log(hazard) 
  equals a linear combination of covariates  + alpha(t) + e where 
  alpha(t)=mu+alpha*log(t) (Weibull random effects)
*/
PROC NLMIXED DATA=jobmult; 
 lambda=EXP(b0+bed*ed+bpres*prestige+bsal*logsal+bseq*seq+e);
 ll=-lambda*duration**(alpha+1)+ event*(LOG(alpha+1)+ 
    alpha*LOG(duration)+LOG(lambda));
 MODEL duration~GENERAL(ll);
 RANDOM e~NORMAL(0,s2) SUBJECT=id;
 PARMS b0=1 bed=0 bpres=0 bsal=0 btime=0 s2=1 alpha=0;
RUN;
/*RANDOM statment says e is normally distributed and subject=id says that e is different 
  random variable for each id variable. 
  PARAMS sets the starting values of the parameters. 
*/


/* The following program is found on page 273 of */
/* "Survival Analysis Using SAS"    */

/*Fix effects models ,
  STRATA is used here so different subgroups may have different baseline hazard functions.*/
PROC PHREG DATA=jobmult NOSUMMARY;
   MODEL duration*event(0)=prestige logsal seq;
   STRATA id;
RUN;
/*Note educ NOT included in this model as fixed effects partial likelihood can only estimate
   coefs for thos covariates that vary across the successive spells for each individual.*/

/*Fix Effects Models is used when: 

1. data do not come from randomized experiment (randome effects preferred in this case)
2. interested in covariates which varies across intervals for each individual
3. most individuals have at least two events 
4. there is a reasonable presumption that the process generating events is invariant over time
5. the unobserved disturbance terim may be correlated with the measured covariates. 

*/


/* The following programs are found on page 276 of */
/* "Survival Analysis Using SAS"    */

/*Method to specify a common origin for all events , a start and stop time must be present */
DATA strtstop;
  SET jobmult;
  RETAIN stop;
  IF seq=1 THEN stop=0;
  start=stop;
  stop=duration+stop;
PROC PRINT;
  VAR id event seq start stop duration;
RUN;

PROC PHREG DATA=strtstop COVSANDWICH(AGGREGATE);
  MODEL (start,stop)*event(0)=prestige logsal ed seq;
  ID id;
RUN;
/*Essentially estimates of counting process model */


/* The following program is found on page 277 of */
/* "Survival Analysis Using SAS"    */

/*Method to specify a common time origin : use only stop time, but essentially stratified:
  or same individual may appear more than once in the same risk set*/
PROC PHREG DATA=strtstop COVSANDWICH(AGGREGATE);
  MODEL stop*event(0)=prestige logsal ed;
  ID id;
  STRATA seq;
RUN;
/*
Disadvantage of using a single time origin: cannot do fixed effects partial likelihood as 
  each individual is treated as a separate stratum in fixed effects partial likelihood. 
*/

/* The following program is found on pages 278-279 of */
/* "Survival Analysis Using SAS"    */

/*Convert duration data into discrete years */
DATA discrete;
  SET jobmult;
  durdis=CEIL(duration);
  DO time=1 TO durdis;
    term=0;
    IF time=durdis THEN term=event;
    OUTPUT;
  END;
RUN;

proc print data=discrete;
run;
/*CEIL Returns the smallest integer that is greater than or equal to the argument,
   fuzzed to avoid unexpected floating-point results.*/

/* The following programs are found on page 279 of */
/* "Survival Analysis Using SAS"    */

/*use SURVEYLOGISTIC and CLUSTER to implement the method of robust standard errors */
PROC SURVEYLOGISTIC DATA=discrete;
  MODEL term(DESC)=prestige logsal ed seq time time*time;
  CLUSTER id;
RUN;
/*
Need to include the time since last event in the covariates
time*time effects need to be specified explicitly while in phreg not */


/* Estimate a proportional hazard model rather than a logit model*/

PROC SURVEYLOGISTIC DATA=discrete;
  MODEL term(DESC)=prestige logsal ed seq time time*time/link=cloglog;
  CLUSTER id;
RUN;


/*Random effects model using GLIMMIX */
PROC GLIMMIX DATA=discrete METHOD=QUAD;
  MODEL term=prestige logsal ed seq time 
    /DIST=BIN SOLUTION LINK=LOGIT;
  RANDOM INTERCEPT / SUBJECT=id;
RUN;
/*METHOD=QUAD forces the procedure to do trum maximum likelihood estimation(used in NLMIXED),
  while the default is pesudo-likelihood estimation(fast but may not accurate for binary outcome).

 Results: the covariance parameter estimates 9.2179 is the estimated variance of e .

*/


/* The following program is found on page 280 of */
/* "Survival Analysis Using SAS"    */

/*Fixed effects realized by LOGISTIC */
PROC LOGISTIC DATA=discrete;
  MODEL term(DESC)=prestige logsal ed seq time time*time;
  STRATA id;
RUN;


/* The following programs are found on page 285 of */
/* "Survival Analysis Using SAS"    */

/*Sensitivity Analysis for the impact informative censoring */

/*Idea: redo the analysis under two extreme assumptions about censored cases

  1. assume censored obs experience events immediately after censored, i.e.,
     censored cases tend to be at high risk of an event 
  2. assume censored obs have longer times to events than anyone else in the sample ,
     i.e., censored cases tend to be at low risk of event */

/*Normal censor */
PROC PHREG DATA=mylib.leaders;	
   MODEL years*lost(0)=manner age start military conflict 
         loginc literacy / TIES=EXACT;
   STRATA region;
RUN;

/*Modify the lost=2 case as years=24, and as a censored case*/
DATA leaders2;
   SET mylib.leaders;
   IF lost=2 THEN years=24;
PROC PHREG DATA=leaders2;	
   MODEL years*lost(0,2)=manner age start military conflict loginc literacy / TIES=EXACT;
   STRATA region;
RUN;
/*In this second case, age effect no longer significant, while age is the only significant 
   covariate for events due to natural death.

We can conclude that treating natural deaths as noninformative censoring has no appreciable 
 effect on the conclusions. 
*/














