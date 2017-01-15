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

-----------------------------------------------------

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

/* The following program is found on page 129 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=myrecid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio;
RUN;



/* The following program is found on page 134 of */
/* "Survival Analysis Using SAS"    */

OPTIONS YEARCUTOFF=1900;
DATA stan;
   INFILE 'c: stan.dat';
   INPUT dob mmddyy9. doa mmddyy9. dot mmddyy9. dls mmddyy9.     
         dead surg m1 m2 m3;
   surv1=dls-doa;
   surv2=dls-dot;
   ageaccpt=(doa-dob)/365.25;
   agetrans=(dot-dob)/365.25;
   wait=dot-doa;
   IF dot=. THEN trans=0; ELSE trans=1;
RUN;

-----------------------------------------------------

/* The following program is found on page 135 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=trans surg ageaccpt;
RUN;
-----------------------------------------------------

/* The following program is found on page 136 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan;
   WHERE trans=1;
   MODEL surv2*dead(0)=surg m1 m2 m3 agetrans wait dot;
RUN;

-----------------------------------------------------

/* The following program is found on page 145 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio 
         / TIES=EXACT;
RUN;

-----------------------------------------------------

/* The following program is found on page 155 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg ageaccpt / TIES=EFRON;
   IF wait>=surv1 OR wait=. THEN plant=0; ELSE plant=1;
RUN;
-----------------------------------------------------

/* The following program is found on page 156 of */
/* "Survival Analysis Using SAS"    */

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
-----------------------------------------------------

/* The following program is found on page 158 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stanlong;
   MODEL (start,stop)*dead2(0)=plant surg ageaccpt /TIES=EFRON;
RUN;

-----------------------------------------------------

/* The following program is found on page 160 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg age / TIES=EFRON;
   IF wait>surv1 OR wait=. THEN plant=0; ELSE plant=1;
   age=ageaccpt+surv1;
RUN;

-----------------------------------------------------

/* The following program is found on page 161 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan;
   MODEL surv1*dead(0)=plant surg logage / TIES=EFRON;
   IF wait>surv1 OR wait=. THEN plant=0; ELSE plant=1;
   logage=LOG(ageaccpt+surv1);
RUN;

-----------------------------------------------------

/* The following programs are found on page 162 of */
/* "Survival Analysis Using SAS"    */

DATA recid;
   INFILE 'c:\recid.dat';
   INPUT week arrest fin age race wexp mar paro prio educ emp1-emp52;
RUN;

PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   DO i=1 TO 52;
      IF week=i THEN employed=emp[i];
   END;
RUN;

-----------------------------------------------------

/* The following program is found on page 163 of */


PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employed=emp[week];
RUN;

-----------------------------------------------------

/* The following programs are found on page 164 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG;
   WHERE week>1;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employed 
         / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employed=emp[week-1];
RUN;


PROC PHREG;
   WHERE week>2;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employ1 
         employ2 / TIES=EFRON;
   ARRAY emp(*) emp1-emp52;
   employ1=emp[week-1];
   employ2=emp[week-2];
RUN;

-----------------------------------------------------

/* The following programs are found on page 165 of */
/* "Survival Analysis Using SAS"    */

DATA recidcum;
   SET recid;
   ARRAY emp(*) emp1-emp52;
   ARRAY cum(*) cum1-cum52;
   cum1=emp1;
   DO i=2 TO 52;
      cum(i)=(cum(i-1)*(i-1) + emp(i))/i;
   END;
PROC PHREG DATA=recidcum;
   WHERE week>1;
   MODEL week*arrest(0)=fin age race wexp mar paro prio employ 
         / TIES=EFRON;
   ARRAY cumemp(*) cum1-cum52;
   EMPLOY=cumemp[week-1];
RUN;


DATA recidcount;
  ARRAY emp(*) emp1-emp52;
  SET recid;
  arrest2=0;
  DO stop=1 TO week;
    start=stop-1;
    IF stop=week THEN arrest2=arrest;
    employed=emp(stop);
    IF week>=1 THEN emplag=emp(stop-1); ELSE emplag=.;
    OUTPUT;
  END;
RUN;

-----------------------------------------------------

/* The following program is found on page 166 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=recidcount;
MODEL (start,stop)*arrest2(0)=fin age race wexp mar paro prio employed /TIES=EFRON;
RUN;

-----------------------------------------------------

/* The following program is found on page 167 of */
/* "Survival Analysis Using SAS"    */

DATA blood;
   INFILE 'c:\blood.dat';
   INPUT deathday status alb1-alb12;
PROC PHREG;
   MODEL deathday*status(0)=albumin;
   ARRAY alb(*) alb1-alb12;
   deathmon=CEIL(deathday/30.4);
   albumin=alb[deathmon];
RUN;

-----------------------------------------------------

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

-----------------------------------------------------

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

-----------------------------------------------------

/* The following program is found on pages 170-171 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=alco;
   MODEL surv*dead(0)=pt;
   time1=0;
   ARRAY time(*) time1-time10;
   ARRAY p(*) pt1-pt10;
   DO j=1 TO 10;
   IF surv > time[j] AND time[j] NE . THEN pt=p[j];
   END;
RUN;

-----------------------------------------------------

/* The following program is found on page 171 of */
/* "Survival Analysis Using SAS"    */

DATA alcocount;
 SET alco;
 time1=0;
 time11=.;
 ARRAY t(*) time1-time11;
 ARRAY p(*) pt1-pt10;
 dead2=0;
 DO j=1 TO 10 WHILE (t(j) NE .);
   start=t(j);
   pt=p(j);
   stop=t(j+1);
   IF t(j+1)=. THEN DO;
      stop=surv;
	dead2=dead;
   END;
   OUTPUT;
END;
PROC PHREG DATA=alcocount;
  MODEL (start,stop)*dead2(0)=pt;
RUN;

-----------------------------------------------------

/* The following program is found on page 173 of */
/* "Survival Analysis Using SAS"    */

ODS GRAPHICS ON;
PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio
         / TIES=EFRON;
   ASSESS PH / RESAMPLE;
RUN;
ODS GRAPHICS OFF;

-----------------------------------------------------

/* The following program is found on page 173 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=recid;
  MODEL week*arrest(0)=fin age prio race wexp mar paro prio /
    TIES=EFRON;
OUTPUT OUT=b RESSCH=schfin schage schprio schrace schwexp schmar 
  schparo schprio;
RUN;

DATA c;
  SET b;
  lweek=log(week);
  week2=week**2;
PROC CORR;
VAR week lweek week2 schfin schage schprio schrace schwexp schmar  
  schparo schprio;
RUN;

-----------------------------------------------------

/* The following program is found on page 178 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio ageweek
         / TIES=EFRON;
   ageweek=age*week;
RUN;
 
-----------------------------------------------------

/* The following program is found on page 180 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=myel;
   MODEL dur*stat(0)=treat / TIES=EXACT;
   STRATA renal;
RUN;

-----------------------------------------------------

/* The following programs are found on page 184 of */
/* "Survival Analysis Using SAS"    */

DATA stan2;
   set stan;
   agels=(dls-dob)/365.25;
RUN;


PROC PHREG DATA=stan2;
   MODEL (ageaccpt,agels)*dead(0)=surg ageaccpt / TIES=EFRON;
RUN;

-----------------------------------------------------

/* The following program is found on page 185 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=stan2;
   MODEL (ageaccpt,agels)*dead(0)=plant surg ageaccpt / TIES=EFRON;
   IF agetrans>=agels OR agetrans=. THEN plant=0;
   ELSE plant=1;
RUN;

-----------------------------------------------------

/* The following program is found on page 186 of */
/* "Survival Analysis Using SAS"    */

ODS GRAPHICS ON;
PROC PHREG DATA=recid PLOTS=S;
   MODEL week*arrest(0)=fin age prio 
         / TIES=EFRON;
   BASELINE OUT=a SURVIVAL=s LOWER=lcl UPPER=ucl;
RUN;
ODS GRAPHICS OFF;

-----------------------------------------------------

/* The following program is found on page 189 of */
/* "Survival Analysis Using SAS"    */

ODS GRAPHICS ON;
PROC PHREG DATA=recid PLOTS(OVERLAY=ROW)=S;
   MODEL week*arrest(0)=fin age prio 
         / TIES=EFRON;
   STRATA fin;
RUN;
ODS GRAPHICS OFF;

-----------------------------------------------------

/* The following program is found on pages 190-191 of */
/* "Survival Analysis Using SAS"    */

DATA covals;
   INPUT fin age prio;
   DATALINES;
0 40 3
;
PROC PHREG DATA=recid;
   MODEL week*arrest(0)=fin age prio / TIES=EFRON;
   BASELINE OUT=a COVARIATES=covals SURVIVAL=s LOWER=lcl   
         UPPER=ucl;
PROC PRINT DATA=a;
RUN;

-----------------------------------------------------

/* The following program is found on page 193 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=recid;
   CLASS educ;
   MODEL week*arrest(0)=fin age prio educ/ TIES=EFRON;
   TEST educ3=educ5;
   CONTRAST 'ed3 vs. ed5' educ 0 1 0 -1 ;
RUN;

CLASS educ(REF='2');

-----------------------------------------------------

/* The following program is found on page 195 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=my.recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio;
   HAZARDRATIO age / UNITS=5 10;
RUN;

-----------------------------------------------------

/* The following programs are found on page 196 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=my.recid;
   MODEL week*arrest(0)=fin age race wexp mar paro prio fin*age;
RUN;


HAZARDRATIO fin / at (age=20 25 30 35 40) CL=PL;

-----------------------------------------------------

/* The following program is found on page 197 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=my.recid;
  MODEL week*arrest(0)=fin age race wexp mar paro prio;
  BAYES;
RUN;



-----------------------------------------------------

---------
Chapter 6
---------

/* The following program is found on page 209 of */
/* "Survival Analysis Using SAS"    */

DATA const;
  SET leaders;
  event=(lost=1);
  type=1;
DATA nat;
  SET leaders;
  event=(lost=2);
  type=2;
DATA noncon;
  SET leaders;
  event=(lost=3);
  type=3;
DATA combine;
  SET const nat noncon;
PROC LIFETEST DATA=COMBINE PLOTS=LLS;
  TIME years*event(0);
  STRATA type;
RUN;

-----------------------------------------------------
/* The following program is found on page 210 of */
/* "Survival Analysis Using SAS"    */

ODS GRAPHICS ON;
PROC LIFETEST DATA=combine PLOTS=H(BW=10);
  TIME years*event(0);
  STRATA type;
RUN;
ODS GRAPHICS OFF;

-----------------------------------------------------
/* The following program is found on page 212 of */
/* "Survival Analysis Using SAS"    */

PROC LOGISTIC DATA=leaders;
   WHERE lost NE 0;
   MODEL lost=years / LINK=GLOGIT;
RUN;

-----------------------------------------------------
/* The following program is found on page 213 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=leaders;
   MODEL years*lost(0)=manner start military age conflict 	
         loginc growth pop land literacy;
   STRATA region;
PROC PHREG DATA=leaders;
   MODEL years*lost(0,1,2)=manner start military age conflict 	
         loginc growth pop land literacy;
   STRATA region;
PROC PHREG DATA=leaders;
   MODEL years*lost(0,1,3)=manner start military age conflict 
         loginc growth pop land literacy;
   STRATA region;
PROC PHREG DATA=leaders;
   MODEL years*lost(0,2,3)=manner start military age conflict 	
         loginc growth pop land literacy;
   STRATA REGION;
RUN;

-----------------------------------------------------
/* The following program is found on page 218 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=combine;
   MODEL years*event(0)=manner start military age conflict 	
         loginc growth pop land literacy / TIES=EFRON;
   STRATA region type;
RUN;

-----------------------------------------------------
/* The following programs are found on page 221 of */
/* "Survival Analysis Using SAS"    */

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


PROC LIFEREG DATA=leaders2; 
  CLASS region;
   MODEL (lower,upper)= manner start military age conflict 	
         loginc literacy region / D=EXPONENTIAL;
RUN;    

-----------------------------------------------------
/* The following program is found on page 225 of */
/* "Survival Analysis Using SAS"    */

DATA leaders3;
   SET leaders;
   lyears=LOG(years+.5);
PROC LOGISTIC DATA=leaders3;
   WHERE lost=1 OR lost=3;
   CLASS region / PARAM=GLM;
   MODEL lost=lyears manner age start military conflict loginc
         literacy region;
RUN;

-----------------------------------------------------
/* The following program is found on page 228 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=leaders;	
   MODEL years*lost(0,2)=manner age start military conflict   
      loginc literacy / TIES=EXACT;
   STRATA region;
RUN;

-----------------------------------------------------
/* The following program is found on page 230 of */
/* "Survival Analysis Using SAS"    */

%CUMINCID(DATA=leaders,
           TIME=years,
           STATUS=lost,
           EVENT=3,
           COMPETE=1 2,
           CENSORED=0)




-----------------------------------------------------

---------
Chapter 7
---------
-----------------------------------------------------
/* The following program is found on page 237 of */
/* "Survival Analysis Using SAS"    */

DATA jobyrs;
   SET jobdur;
   DO year=1 TO dur;
      IF year=dur AND event=1 THEN quit=1; 
      ELSE quit=0;
      OUTPUT;
   END;
RUN;

-----------------------------------------------------
/* The following program is found on page 238 of */
/* "Survival Analysis Using SAS"    */
 
PROC LOGISTIC DATA=jobyrs;
  CLASS year / PARAM=GLM;
  MODEL quit(DESC)=ed prestige salary year;
RUN;


-----------------------------------------------------
/* The following program is found on page 244 of */
/* "Survival Analysis Using SAS"    */
 

DATA rankyrs;
   INFILE 'c:rank.dat';
   INPUT dur event undgrad phdmed phdprest art1-art10 
         cit1-cit10 prest1 prest2 jobtime;
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


-----------------------------------------------------
/* The following program is found on page 245 of */
/* "Survival Analysis Using SAS"    */


PROC LOGISTIC DATA=rankyrs;
   MODEL promo(DESC)=undgrad phdmed phdprest articles citation 
         prestige year year*year;
RUN;


-----------------------------------------------------
/* The following programs are found on page 251 of */
/* "Survival Analysis Using SAS"    */

DATA jobyrs2;
   SET jobdur;
   DO year=1 TO dur;
     IF year=dur THEN outcome=event;
     ELSE outcome=0;
     OUTPUT;
   END;
RUN;

PROC LOGISTIC DATA=jobyrs2;
  MODEL outcome(REF='0')=ed prestige salary year 
    / LINK=GLOGIT;
RUN;

-----------------------------------------------------
/* The following programs are found on page 253 of */
/* "Survival Analysis Using SAS"    */


PROC LOGISTIC DATA=jobyrs2;
   WHERE outcome NE 2;
   MODEL outcome(DESC)=ed prestige salary year;
RUN;

PROC LOGISTIC DATA=jobyrs2;
   WHERE outcome NE 1;
   MODEL outcome(DESC)=ed prestige salary year;
RUN;

-----------------------------------------------------
/* The following program is found on page 262 of */
/* "Survival Analysis Using SAS"    */

PROC SORT DATA=jobmult;
   BY seq;
PROC PHREG DATA=jobmult;
   MODEL duration*event(0)=prestige logsal ed;
   BY seq;
RUN;

-----------------------------------------------------
/* The following program is found on page 265 of */
/* "Survival Analysis Using SAS"    */

PROC SORT DATA=jobmult;
  BY id seq;
DATA joblag;
  SET jobmult;
  durlag=LAG1(duration);
PROC PHREG DATA=joblag;
  WHERE seq = 2;
  MODEL duration*event(0)=prestige logsal ed durlag;
RUN;

-----------------------------------------------------
/* The following program is found on page 267 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=jobmult COVSANDWICH(AGGREGATE);
   MODEL duration*event(0)=prestige logsal ed seq;
   ID id;
RUN;

-----------------------------------------------------
/* The following programs are found on page 269 of */
/* "Survival Analysis Using SAS"    */

DATA job;
  SET my.jobmult;
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

-----------------------------------------------------
/* The following program is found on page 271 of */
/* "Survival Analysis Using SAS"    */

PROC NLMIXED DATA=jobmult; 
 lambda=EXP(b0+bed*ed+bpres*prestige+bsal*logsal+bseq*seq+e);
 ll=-lambda*duration**(alpha+1)+ event*(LOG(alpha+1)+ 
    alpha*LOG(duration)+LOG(lambda));
 MODEL duration~GENERAL(ll);
 RANDOM e~NORMAL(0,s2) SUBJECT=id;
 PARMS b0=1 bed=0 bpres=0 bsal=0 btime=0 s2=1 alpha=0;
RUN;

-----------------------------------------------------
/* The following program is found on page 273 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=jobmult NOSUMMARY;
   MODEL duration*event(0)=prestige logsal seq;
   STRATA id;
RUN;

-----------------------------------------------------
/* The following programs are found on page 276 of */
/* "Survival Analysis Using SAS"    */

DATA strtstop;
  SET jobmult;
  RETAIN stop;
  IF seq=1 THEN stop=0;
  start=stop;
  stop=duration+stop;
PROC PRINT;
  VAR id event seq start stop;
RUN;

PROC PHREG DATA=strtstop COVSANDWICH(AGGREGATE);
  MODEL (start,stop)*event(0)=prestige logsal ed seq;
  ID id;
RUN;

-----------------------------------------------------
/* The following program is found on page 277 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=strtstop COVSANDWICH(AGGREGATE);
  MODEL stop*event(0)=prestige logsal ed;
  ID id;
  STRATA seq;
RUN;

-----------------------------------------------------
/* The following program is found on pages 278-279 of */
/* "Survival Analysis Using SAS"    */

DATA discrete;
  SET jobmult;
  durdis=CEIL(duration);
  DO time=1 TO durdis;
    term=0;
    IF time=durdis THEN term=event;
    OUTPUT;
  END;
RUN;

-----------------------------------------------------
/* The following programs are found on page 279 of */
/* "Survival Analysis Using SAS"    */

PROC SURVEYLOGISTIC DATA=discrete;
  MODEL term(DESC)=prestige logsal ed seq time time*time;
  CLUSTER id;
RUN;

PROC GLIMMIX DATA=discrete METHOD=QUAD;
  MODEL term=prestige logsal ed seq time 
    /DIST=BIN SOLUTION LINK=LOGIT;
  RANDOM INTERCEPT / SUBJECT=id;
RUN;

-----------------------------------------------------
/* The following program is found on page 280 of */
/* "Survival Analysis Using SAS"    */

PROC LOGISTIC DATA=discrete;
  MODEL term(DESC)=prestige logsal ed seq time time*time;
  STRATA id;
RUN;

-----------------------------------------------------
/* The following programs are found on page 285 of */
/* "Survival Analysis Using SAS"    */

PROC PHREG DATA=leaders;	
   MODEL years*lost(0)=manner age start military conflict 
         loginc literacy / TIES=EXACT;
   STRATA region;
RUN;

DATA leaders2;
   SET leaders;
   IF lost=2 THEN years=24;
PROC PHREG DATA=leaders2;	
   MODEL years*lost(0,2)=manner age start military conflict 				   loginc literacy / TIES=EXACT;
   STRATA region;
RUN;















