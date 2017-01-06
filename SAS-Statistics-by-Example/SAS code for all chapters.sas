

/*Chapter 1 An introduction to SAS */
/*Print the imported data out*/
proc print data=SampleData;
run;

/*Imoport txt files */
data Sample2;
infile 'C:\Users\Qianqian Shan\Dropbox\SAS Statistics by example\delim.txt';
length Gender $ 1;
input ID Age Gender $;
run;


/*Add a title that prints accross the top of every page of output*/
title 'Listing of Data Set Sample2';
proc print data=Sample2;
run;


/*Import csv files
  the dsd at the end of infile:
  specify that two consecutive commas 
  represent a missing value and 
  the default delimiter is a comma*/
data Sample3;
infile 'C:\Users\Qianqian Shan\Dropbox\SAS Statistics by example\comma.csv' dsd;
length Gender $ 1;
input ID Age Gender $;
run;

title 'Listing of Data Set Sample3';
proc print data=Sample3;
run;


/*Chapter 2 Descriptive Statistics-continuous variables */

/*n:number of nonmissing obs
  nmiss: number of obs with missing values
  maxdec:maximum number of decimal places
  clm:95% confidence limit on the mean
  cv: coefficient of variation */
title Descriptive Statistics for SBP and DBP;
proc means data=example.Blood_Pressure n nmiss mean std median maxdec=3;
var SBP DBP;
run;

/*use class to see statistics for each level*/
title Descriptive Statistics for SBP and DBP by drug;
proc means data=example.Blood_Pressure n nmiss mean std median maxdec=3;
class drug;
var SBP DBP;
run;

/*use printalltypes to print both the grand mean and the mean by class*/
title Descriptive Statistics for SBP and DBP:print all types;
proc means data=example.Blood_Pressure n nmiss mean std median printalltypes maxdec=3;
class drug;
var SBP DBP;
run;

/*compute the confidence interval using clm and standard error
to decide how well the sample mean estimates the population mean*/
title Descriptive Statistics for SBP and DBP: 95% confidence interval;
proc means data=example.Blood_Pressure n nmiss mean std median clm stderr printalltypes maxdec=3;
class drug;
var SBP DBP;
run;

/*simple print the data out for reference*/
proc print data=example.Blood_Pressure;
run;

/*Use univariate procedure to produce histogram and probability plots
ID statement not necessary but useful to specify a variable that 
identify each observation*/
title 'Produce histograms and probability plots using proc univarate';
proc univariate data=example.Blood_Pressure;
id Subj;
var SBP DBP;
histogram /midpoints=70 to 100 by 5 normal;
probplot / normal(mu=est sigma=est);
run;

/*sgplot procedure is used to produce histogram*/
title Using SGPLOT to produce histogram;
proc sgplot data=example.Blood_Pressure;
histogram SBP;
run;

title Using SGPLOT to produce boxplot;
proc sgplot data=example.Blood_Pressure;
hbox SBP;
/*or vbox to produce vertical box plot*/
run;


title 'Generate new data with outliers';
data Blood_Pressure_Out;
set example.Blood_Pressure(keep=Subj SBP);
if Subj=5 then SBP=200;
else if Subj=55 then SBP=180;
run;

/*produce box plot and label the outliers*/
title 'Display outliers in box plot';
proc sgplot data=Blood_Pressure_Out;
hbox SBP /datalabel=SBP;
run;

title 'Display multiple box plots for each category';
proc sgplot data=example.Blood_Pressure;
hbox SBP/category=Drug;
run;

/*SGPLOT can be used to generate histograms,box plots,scatter plots
SGSCATTER displays several plots on a single page*/





/*Chapter 3 Descriptive Statistics-Categorical variables */

title 'computing frequencies and percentages using proc freq';
proc print data=example.Blood_Pressure;
run;

proc freq data=example.Blood_Pressure;
tables Gender Drug;
run;

/*demonstrate the noncum table option*/
proc freq data=example.Blood_Pressure;
tables Gender Drug /nocum missing;
/*the missing option tells proc to treat missing values as a valid category
and include them in the body of the table, not necessary */
run;

/*compute frequency on continuous variables
this one will display each unique variable value
in one row */
proc freq data=example.Blood_Pressure;
tables SBP /nocum missing;
run;

/*can use formats to group the continuous variables:
create a format to associate any text to one or more values*/
proc format;
value $gender 'M'='Male'
              'F'='Female';
value sbpgroup low-140='Normal'
               141-high='High';
value dbpgroup low-80='Normal'
               81-high='High';
run;

/* use the format defined above the group continuous data
the period after $gender is required for a statement
format for character variables must start with $
format for numeric variables must not start with $
format names cannot end with numbers */
proc freq data=example.Blood_Pressure;
tables Gender SBP DBP/nocum;
format Gender $gender.
       SBP    sbpgroup.
	   DBP    dbpgroup.;
run;


/*generate a bar chart using sgplot*/
proc print data=store;
run;
proc sgplot data=store;
vbar Region;
run;

/*use ods to create pdf/html/rtf files*/
ods listing close; /*to avoid getting two outputs*/
ods pdf file='c:\Users\Qianqian Shan\Dropbox\SAS Statistics by example\bar.pdf' style=journal;
/*could also use style statistical*/
title 'Generate a bar chart using proc sgplot';
proc sgplot data=store;
vbar region;
run;
quit; /*quit the ods*/
ods pdf close;
ods listing;

/*generate html file*/
ods listing close; /*to avoid getting two outputs*/
ods html body='c:\Users\Qianqian Shan\Dropbox\SAS Statistics by example\bar.htm' style=HTMLBlue;
title 'Generate a bar chart using proc sgplot';
proc sgplot data=store;
vbar region;
run;
quit; /*quit the ods*/
ods html close;
ods listing;


/*create cross-tablation table by proc freq*/
/*tje row and colums order are determined by the internal value*/
title 'Demonstrate cross tabulation';
proc freq data=store;
tables Gender*Region; /*tables rowname * colname */
run;

/*change the order of values*/
proc format;
value $region 'North'='1 North'
              'East'='2 East'
			  'South'='3 South'
			  'West'='4 West';
run;

title 'Demonstrate cross tabulation with ordered names';
/*need to tell sas the table is ordered by the formatted values 
by order=formatted*/
proc freq data=store order=formatted;
tables Gender*Region; /*tables rowname * colname */
format Region $region.;
run;





/*Chapter 4: Descriptive statistics-bivariate associations*/

title 'Creating a scatter plot using proc gplot';
symbol value=dot; /*a global graphics instruction to tell proc to use dot*/
proc gplot data=store;
plot Book_Sales*Music_Sales;  /*plot yvariable * xvariable*/
run;
quit;

/*Scatter plot to show different genders in different colors*/
goptions reset=all; /*graphics options stay in effect until you change them so it's good to reset it*/
title 'Creating a scatter plot sorted by Gender';
title2 'Gender is used to sort the dots';
symbol1 color=red value=diamond;/*a global graphics instruction to tell proc to use dot*/
symbol2 color=blue value=circle; /*or square,hash,circle,dot...*/
proc gplot data=store;
plot Book_Sales*Music_Sales=Gender;  /*plot yvariable * xvariable*/
run;
quit;


/*Using SGPLOT to produce scatter plot*/

title 'Using sgplot to produce a scatter plot';
proc sgplot data=store;
scatter x=Book_Sales y=Music_Sales /group=Gender makeattrs=(symbol=trianglefilled);
run;
quit;

/*Create multiple scatter plots on a single page using proc sgscatter*/
title 'Demonstrating the plot statement by sgscatter';
proc sgscatter data=store;
plot Book_Sales*Music_Sales Total_Sales*Electronics_Sales;
/*Just list the multiple plots by y*x y1*x1 ...*/
run;


title 'Demonstrate the compare statement by sgscatter';
proc sgscatter data=store;
compare y=Total_Sales x=(Book_Sales Music_Sales Electronics_Sales)/group=Gender;
run;

title 'Demonstrate the matrix statement by sgscatter';
proc sgscatter data=store;
matrix Book_Sales Music_Sales Electronics_Sales/diagonal=(histogram);
/*diagonal option is used to contain histogram on the diagonal, could also be normal,kernel etc*/
run;


/*Chapter 5: Inferential statistics one-sample tests*/
/*Proc TTEST can be used to perform t-test.two sample t test unpaired and paired*/

/*One sample t-test using PROC TTEST (could also use UNIVARIATE, later)*/
proc ttest data=exercise h0=50 sides=2 alpha=0.5;
var Age;
run;

/*One sample t-test using PROC UNIVARIATE*/
title 'Conducting one sample t-test with a given population mean';
proc univariate data=exercise mu0=50; /*Note it's mu0 instead of mu*/
var Age;
id Subj;
run;

/*Test for normality*/
/*Not always useful(for large sample sizes)*/
title 'Testing if a variable is Normally distributed or not';
proc univariate data=exercise normal;
/*Can also specify other distributions to test: normal(by default)
beta,exponential,gamma,lognormal,two parameter Weibull,three parameter weibull...*/
var Age;
probplot /normal (mu=est sigma=est);
run;

/*Chapter 6: Inferential test: two sample tests*/

/*Two Sample ttest by ttest, under the assumption of equal/unequal variances*/
title 'Conduct two sample t-test by TTEST';
proc ttest data=example.blood_pressure;
class Gender;
var SBP DBP;
run;

/*Use ods graphics to produce all plots for the two sample t test*/
ods graphics on;
proc ttest data=example.blood_pressure
plots(unpack shownull)=all; 
/*use unpack to show all plots
use shownull to display a vertical reference line at the null hypothesis value
(default is zero)
and show the 95% CI by default*/
class Gender;
var SBP DBP;
run;
ods graphics off;


/*Conduct a paired t-test for paired data */
title 'Demonstrate a Paired T-test';
proc ttest data=reading;
paired After*Before;
run;
 


/*Chapter 7 Inferential statistics- compare more than two means-ANOVA*/
/*GLM procedure is used */
ods graphics on;
title 'Running a One-Way ANOVA using PROC GLM';
proc glm data=store plots=diagnostics; 
/*use plots(unpack)=diagnostics to show
diagonostics plots seperately, diagnostics can show an additional panel of plots*/
class Region;
model Electronics_Sales=Region/ss3; /*Show only the type III SS*/
means Region /hovtest welch; /*the class model means order matters*/
run;
quit;
ods graphics off;

/*Note: ANOVA is robust with regard to the assumption of homogeneity of variance,
especially if you have balanced design*/

/*If you decide there is a violation of the homo variance assumption,can consisder
non parametric alternatives to ANOVA, or use WELCH option on the mean statement: 
gives weights to the cell means based on the inverse of the variance*/ 

/*Once you ar satisfied with the test assumptions, you can continue to consider the ANOVA results
and do multiple comparision tests */

/*MULTIPLE COMPARISON: USE TUKEY STUDENT NEWMAN KEULS (SNK) OPTION FOR MEAN*/
/*MULTIPLE COMPARISON: USE PDIFF OPTION FOR LSMEANS*/

title 'Request multiple comparison tests';
proc glm data=store;
class Region;
model Electronics_Sales=Region/ss3;
means Region/snk;
lsmeans Region /pdiff adjust=tukey;
/*Other adjust method includes:BON,DUNNETT(DUNNETT'S TEST)...,the default is TUKEY*/
run;
quit;




/*TWO WAY FACTORIAL DESIGNs*/
/*two classes in the class statement
model statement has one more factor after the = sign*/
title 'Perfom two-way factorial design';
proc glm data=store;
class Region Gender;
model Electronics_Sales=Region |Gender /ss3;
/*can also use : model Electronics_Sales=Region Gender Region*Gender; 
use @2 to place a limit on the interaction terms if necessary*/
lsmeans Region|Gender;
run;
quit;


/*Another example with significant interaction*/
title 'Perfom two-way factorial design with significant interaction';
proc glm data=store;
class Region Gender;
model Music_Sales=Region |Gender /ss3;
/*can also use : model Electronics_Sales=Region Gender Region*Gender; 
use @2 to place a limit on the interaction terms if necessary*/
lsmeans Region|Gender;
lsmeans Region*Gender/slice=Region; 
/*compare between specific regions and genders for the significant interaction*/
lsmeans Region*Gender/slice=Gender; 
run;
quit;


/*RANDOMIZED BLOCK DESIGN*/
/*It's similar to the previous factorial design*/
proc glm data=store;
class Region Gender;
model Book_Sales=Gender Region /ss3;
lsmeans Gender;
run;
quit;
/*In this model, Region serves as a block and assume no interaction between block and factor*/


/*Chapter 8 Correlation and Regression */

/*Produce Pearson Correlations*/
ods graphics on;
title 'Compute Pearson Correlation Coef';
/*nosimple tells the procedure that you don't want the defalt output of means and sd
for ach of the variables in the VAR and WITH
RANK-order the correlations from largest tot smallest by abs values*/
proc corr data=exercise nosimple rank plots=matrix(nvar=all); /*or nvar=n*/
var Rest_Pulse Max_Pulse Run_Pulse Age;
with Pushups;
run;
ods graphics off;


/*Or use plots option to see each scatter plot seperately
usually up to 10 plots*/
ods graphics on;
title 'Compute Pearson Correlation Coef';
/*nosimple tells the procedure that you don't want the defalt output of means and sd
for ach of the variables in the VAR and WITH
RANK-order the correlations from largest tot smallest by abs values*/
proc corr data=exercise nosimple rank plots(only)=scatter(nvar=all ellipse=none);
/* only says you only want separate bivariate pltos for each variable pair
or nvar=n , sas draw 95% prediction ellipse on the scatter plots by default, 
can use ellipse =none to eliminate it*/
var Rest_Pulse Max_Pulse Run_Pulse Age;
with Pushups;
run;
ods graphics off;

/*Or use SGSCATTER to produce as many plots as you want*/




/*Generate Correlation Matrix*/
ods graphics on;
title 'Compute Pearson Correlation Matrix';
/*omit the WITH statement*/
proc corr data=exercise nosimple plots=matrix(histogram);
var Pushups Rest_Pulse Max_Pulse Run_Pulse Age;
run;
ods graphics off;


/*Create HTML ouput with data tips*/
ods graphics on/imagemap=on;
ods listing close;
ods html gpath='C:\Users\Qianqian Shan\Dropbox\SAS Statistics by example'
         path='C:\Users\Qianqian Shan\Dropbox\SAS Statistics by example'
		 file='scatter.html' 
		 style=science;/*journal or statistical,science*/
		 /*save html fine in with the file name and in the path 
		 save each plot separately in gpath with names given automatically */
title 'Save Person Correlation Coef in a html file';
proc corr data=exercise nosimple plots(only)=scatter(ellipse=none);
var Rest_Pulse Max_Pulse Run_Pulse Age;
with Pushups;
run;
ods html close;
ods graphics off;
ods listing;

/*Generate Spearman nonparametric correlations*/
title 'Compute Spearman Rank corrlation and Pearson correlation';
proc corr data=exercise nosimple spearman pearson; /*specify both types here*/
var Rest_Pulse Max_Pulse Run_Pulse Age;
with Pushups;
run;

/*Run a simple linear regression model*/
title 'A simple linear regression model';
proc reg data=exercise plots(unpack)=diagnostics;
model Pushups=Rest_Pulse;
run;
quit;

/*use ODS statistical graphics to investigate influential obs:
Run the regression model with and without a given obs and see the effect
on the predicted value or on the betas */

ods listing close;
title 'Display Influential Observations';
proc reg data=exercise plots(only)=(cooksd(label),rstudentbypredicted(label));
id Subj;
model Pushups=Rest_Pulse /influence r;
run;
quit;
ods listing;
/*INFLUENCE option gives the statistics taht show how much each obs changes
 aspects of the regression depending on whether that observation is included

R option gives more details about the residuals as well as the value of Cook's D statistic

Cook's D measures the effect of each data point on the predicted value 

Label option in PLOTS places values of the id on the plot to identify data points 
that exceed the standard cutoff.
*/


/*Meothod One: Use regression Equation to do prediction */
/*Method: Create a data set that contains the values of the covariates*/
title 'Create data for prediction and combined it with the original data';
data need_prediction;
input Rest_Pulse @@;
datalines;
50 60 70 80 90
; /*semicolon has to be in another line*/
data combined;
set exercise need_prediction;
run;
proc print data=combined;
run;

title 'use proc reg to compute predicted values';
proc reg data=combined;
model Pushups=Rest_Pulse /p;
id Rest_Pulse;
run;
quit;

/*p option in the model statement prints the predicted values for each obs

id ads the varaiable Rest_Pulse to the output so we can know what x the predicted value
corresponds to */

/*Method two of making prediction write the model parameters into a data set and use them to compute
predicted values by PROC SCORE
*/

title 'Save the parameters';
proc reg data=exercise noprint outest=betas;
model Pushups=Rest_Pulse;
run;
quit;
proc print data=betas noobs;
run;
/*NOPRINT tells the procedure that no need for any printed output

OUTEST=(output estimates) tells that you want the regression parameters to be
written to a data set called betas

NOOBS tells that don't print the obs column, which is like an index
*/

title 'use PROC SCORE to make prediction';
proc score data=need_prediction score=betas out=predictions type=parms;
var Rest_Pulse;
run;
proc print data=predictions noobs;
run;
/*score= tells the procedure the name of the data sets that contains the regression param

out= tells the name of the data set that will contain the predicted values

type= tells the type of data it's dealing with, can also be factors...*/






/*Chapter 9 Multiple Regression */

/*Run a multiple regression model*/
title 'Run a multiple regression model';
proc reg data=exercise;
model Pushups=Age Max_Pulse;
run;
quit;


/*Model selection: RSQUARE
the process runs 2^n-1 possible models for n predictor variables*/

title 'Demonstrate the RSQUARE selection';
proc reg data=exercise;
model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=rsquare cp adjrsq best=2;
run;
quit; 
/* SELECTION indicates rsquare method is used

CP indicates Mallows' Cp, which can be used to decide if you have too few 
or too many predictor variables

ADJRSQ is adjusted R square values

BEST= shows the number of models to display for each case */


title 'Demonstrate the RSQUARE selection with separate plots';
proc reg data=exercise plots(unpack)=(diagnostics residualplot);
model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=rsquare cp adjrsq best=2;
run;
quit; 

title 'Demonstrate the RSQUARE selection with partial separate plots';
proc reg data=exercise plots=(diagnostics residualplot(unpack));
model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=rsquare cp adjrsq best=2;
run;
quit; 

/*Model Selection: Mallows Cp and Hocking's Criteria*/
/*Mallows suggested to choose the first model in which Cp is less than or equal
the number of parameters. UAed when using a regression to make predictions among the predictors
and denpendent variables. 

Hocking proposed to choose the first model where Cp is less than 2*p-p_full+1, where p_full
is the number of paramters in a full model
*/

title 'Generate plots of Rsquare Adjusted Rsquare and Cp';
proc reg data=exercise plots(only)=(rsquare adjrsq cp);
model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=rsquare cp adjrsq;
run;
quit;
/*plots(only)=(rsquare adjrsq cp) prints plots of rsquare, adjusted rsquare and Cp
   vs. Number of parameters 

Mallows criteria will select more parsimonious models in many cases, i.e.,a simpler model
with fewer parameters */


/*Model Selection: Forward, Backward, Stepwise */
/*
Forward: Start with a single best variable and adds variables one at a time until 
         the p value of the variable being entered is larger than the specified value
Backward: Start with the full model and removes one at a time until all variables being 
          considered for removal have p value smaller than a specified value
Stepwise: Almost the same as Forward, but a variable has been added may be removed later
          if the p value is larger than the specified one later 
Two options: SLENTRY(significance level for entering), SLSTAY(significance level to stay)

Default values: 
Forward: SLENTRY=0.5
Backward:SLSTAY=0.1
Stepwise: SLENTRY=0.15, SLSTAY=0.15

Default values can be changed by add an option in the model statement.
*/
title 'Forward, backward and stepwise methods';
proc reg data=exercise;
Forward: model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=forward;
Backward: model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=backward;
Stepwise: model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/selection=stepwise;
run;
quit;
/*The words before model (Forward:,Backward: ...) serve as label of each model to be 
  included in the output

selection is written as an option for the model

the p value always represent the contribution to each varible after all the other variables 
 have been accounted for */



/*Force selected variables into a model */
/*First put the variable to be included as the first variable after model= ,then add 
  include=1 in option*/
title 'Forcing variables into a stepwise model';
proc reg data=exercise;
model Pushups=Max_Pulse Age Rest_Pulse Run_Pulse/selection=stepwise include=1;
run;
quit;


/*Create dummy variables for regression when including categorical variables*/

data dummy;
set Store;
*Create dummy variables for Gender;
if Gender='Male' then Male=1;
else if Gender='Female' then Male=0;

*Create dummy variables for Region;
if Region not in ('East' 'North' 'South' 'West') then 
call missing(North,East,South);
else if Region='North' then North=1;
else North=0;
if Region='East' then East=1;
else East=0;
if Region='South' then South=1;
else South=0;
run;
proc print data=dummy(obs=10) noobs;
var Region Gender Male North East South;
run;
/*CALL MISSING sets each of the three variables in the bracket to a missing value */

/*Run multiple regression with Dummy variables */
title 'Run Regression with Dummy variables';
proc reg data=Dummy;
model Music_Sales=Total_Sales Male North East South;
run;
quit;


/*Detect Collinearity*/
/*If two or more variables are correlated and they both entered into the regression 
model , you might get very unstable estimates of slopes */
title 'Use VIF to detect collinearity';
proc reg data=exercise;
model Pushups=Age Rest_Pulse Max_Pulse Run_Pulse/VIF;
run;
quit;
/*Variance Inflation Factor is used to assess collinearity: Run a regression of all 
  other predictor variables to predict the one in question and VIF=1/(1-R_i^2)

  VIF greater than 10 will be considered large */

/*Influential Observations in Multiple Regression Model */
title 'Detect Influential Observations in Multiple Regression';
proc reg data=exercise plots(label only)=(cooksd rstudentbypredicted dffits dfbetas);
id Subj;
model Pushups=Age Max_Pulse Run_Pulse /influence;
run;
quit;
/*Four standards are saved here
Cooksd - measure the effect on the predicted value
Rstudentbypredic -externally studentized residuals by predicted values
DFFITS -the difference in the overall effect on the betas
DFBETAS -the difference on each beta 

See https://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/
viewer.htm#statug_reg_sect040.htm
for the formulas to calculate these values
*/





/*Chapter 10 Categorical Data */
/*Compare the proportions using the chi-square test and Fisher's exact test*/
/*Relative Risk and Odds Ratio are commonly used to describe the results of cohort
  and case-control studies */

title 'Compare Proportions';
proc freq data=risk;
tables Gender*Heart_Attack/all;
run;
/*tables row-variablse*col-variablse/ options 
Commmonly Used Options on the Tables statement:
AGREE: Kappa Coefficient of agreement
ALL: A quick way to get CHISQ,MEASURES,and CMH
CHISQ: chisquare,Fisher's exact(for 2*2) and several other tests
CMH: Cochran-Mantel-Haenszel statistics
FISHER: Fisher's exact test
MEASURES: various meansures of association, such as Pearson and Spearman coefficients
RELRISK: Relative risk(risk ratio)
TREND: Cochran-Armitage test for trend 
*/


/*If more than 20% of the cells have expected values that are less than 5, 
  computing chi-square might not be approciate 

For 2*2 tables: use Yate's correction for continuity or Fisher's exact test 
For larger tables: request Fisher's exact test or collapse cells in the table 
 so that the cell frequencies are largers */

data small_counts;
input Group $ Outcome $ @@;
datalines;
A Good A Good A Good A Poor A Good A Good
B Poor B Poor B Good B Poor B Poor B Poor
;
title 'Tables with Small Expected Frequencies';
proc print data=small_counts;
run;
proc freq data=small_counts;
tables Group*Outcome/chisq;
run; 


/*Compute chi-square from Frequency table*/
/*Applies when we have a table data already */

data freq;
input Treatment $ Outcome $ Count;
datalines;
Placebo Sick 30
Placebo Well 10
Drug Sick 15
Drug Well 40
;
title 'Compute Chi-square from Frequency Data';
proc freq data=freq;
tables Treatment*Outcome/chisq;
weight Count;
run;
/*WEIGHT specifies the frequencies for each combination*/



/*Use MACRO */
/***********************************************************
Macro CHISQ
Purpose: To compute chi-square (and any other valid
         PROC FREQ TABLES options) from frequencies in a
         2 x 2 table.
Sample Calling Sequencies;
   %CHISQ(10,20,30,40)
   %CHISQ(10,20,30,40,OPTIONS=CMH)
   %CHISQ(10,20,30,40,OPTIONS=CHISQ CMH)
************************************************************/
%macro chisq(a,b,c,d,options=chisq);
   data chisq;
      array cells[2,2] _temporary_ (&a &b &c &d);
      do row = 1 to 2;
         do Col = 1 to 2;
            Count = cells[Row,Col];
            output;
         end;
      end;
   run;
   proc freq data=chisq;
      tables Row*Col / &options;
      weight Count;
   run;
%mend chisq;
/*The above is the CHISQ macro. 
Another way to call it is write the following command:
%include 'c:\..\.sas */

/*Example to call the macro to do calculation */
title;
%chisq(30,10,15,40,options=all) /*No Semi colon here, in order topleft, topright,bottomleft,bottomright */


/*Compute Coefficient Kappa-A test of Agreement */
data kappa;
input Rater1 : $1. Rater2 : $1. @@;
datalines;
Y Y N N Y N N Y Y Y Y Y Y Y N N N N N N Y Y Y N Y Y N N N Y N N N N N N
Y Y Y Y N N N N Y N Y Y Y Y N N N N N N N N Y Y N N Y Y N N 
;
proc print data=kappa;
run;

title 'Compute Coefficient Kappa';
proc freq data=kappa;
tables Rater1*Rater2/agree;
test kappa;
run;
/*Override this default vlue of the storage length of a character variables:
Method 1: follow the variable with a colon, and followed by an informat
Method 2: include an INFORMAT statement before INPUT statement, which associates
          each variable with the informat you want to use */
/*A Kappa value of 1 indicates perfect agreement, and a value 0 indicates that the 
  agreement is attributable to chance */


/*Compute Tests for Trends */

data trend;
input Outcome $ Dose Count;
datalines;
Success 1 8
Success 2 8
Success 3 10
Success 4 15
Failure 1 12
Failure 2 12
Failure 3 10
Failure 4 5
;
title 'Compute tests for trends';
proc freq data=trend;;
tables Outcome*Dose/cmh trend;
weight Count;
run;
/*Cochran-Mantel-Haenszel chisquare and Cochran-Armitage test for trend are done.

 P values from both tests indicate that success rate increases with increasing doses */


/**/
data bloodtype;
input Type $ Count @@;
datalines;
O 88 A 76 B 24 AB 12
;

title 'Compute chi-square for one-way table';
proc freq data=bloodtype;
tables Type/testp=(.4 .04 .11 .45);
weight Count;
run; 
/*TESTP= options tells the expected proportions either as numbers less than 1 
         (adds up tp 1)or percentages (adds up to 100%)

  use the TESTP option to test if the available bloodtype data has the same disttribution 
  as was found in the whole country */



/*Chapter 11 Binary Logistic Regression */
title 'Logstic Regression with One Categorical Predictor Variable';
proc logistic data=risk;
class Gender (param=ref ref='F');
model Heart_Attack (event='Yes')=Gender/clodss=pl ;
run;
quit;
/* PARAM=REF use reference coding paramterization for dummy variables, 
   ref='F' chooses F as reference level
   Another paramterization is effect coding, which uses the mean level of the 
     predictor variabls as the reference level; 
   In comparison, reference coding uses a single value of the predictor variables
     as a reference as shown in this example

EVENT='Yes' shows that you want to predict the probability of having a heart-attack(Yes)
 if it's not specified, the first value appeared in alphabetical or numerical order will be 
 used .
Note: If event is a numeric variable, single quote is still needed. i.e. EVENT='1'

CLODDS=pl indicates that the confidence limits of odds are based on profile likelihood,
 the default is Wald type interval, this is shown as the last part of the output
*/


/*Run a logistic regression model with one continuous predictor variable*/
title 'Logistic Regression with one continuous variable';
proc logistic data=risk;
model Heart_Attack(event='Yes')=Chol/clodds=pl;
units chol=10; /*unit or units ,both work*/
run;
quit;


/*Use a format to create a categorical variable from a continuous variable */

proc format;
value cholgrp low-200 ='Low to Medium'
              201-high ='High';
run;
title 'Use format to create categorical variable from continuous variable';
proc logistic data=risk;
class Chol (param=ref ref='Low to Medium');
model Heart_Attack(event='Yes') =Chol/clodds=pl;
format Chol cholgrp.;
run;
quit;
/*The reference level in the CLASS statement is based on the formatted variables*/ 


/*Use combination of categorical and continuous variables */
title 'Logistic regression with categorical and continuous variables';
proc print data=risk;
run;
ods graphics on;
proc logistic data=risk;
class Age_Group(ref='1:< 60') Gender(ref='F') /param=ref;
model Heart_Attack(event='Yes')=Gender Age_Group Chol /clodds=pl;
units Chol=10;
run;
quit;
ods graphics off;
/*the name in the single quote of the REF has to be exactly the same as shown in the data

  If there are multiple class variables and we want to use reference coding for all of them,
   we can use param=ref in the options */


/*Logistic regression with interactions */
title 'Run logistic model with interactions';
proc logistic data=risk plots(only)=(roc oddsratio);
class Gender(param=ref ref='F') Age_Group(param=ref ref='1:< 60');
model Heart_Attack(event='Yes')=Gender |Age_Group |Chol @2/selection=backward slstay=0.1 clodds=pl;
units Chol=10;
oddsratio Chol;
oddsratio Gender;
oddsratio Age_Group;
run;
quit;
/*In the backward elimination method, main effects cannot be removed from the model if those 
   effects are involved in an significant interaction (hierarchical model) 

 ROC shows the relationship between a false-positive rate and the sensitivity of the test. 
   False-positive rate: the porportion of time an obs without the outcome is predicted by 
      the model to have the outcome. 
   Sensitivity: the proportion of obs that have the outcome that is predicted by the model 
      to have the outcome. 
 Ideally, models with high sensitivity and low false-positive rate is preferred. 
 The more area under the ROC curve, the more ideal the model. 

When interactions involved in the model, the LOGISTIC procedure will not compute oddsratio 
  automatically, so we need to speicify it explicitly. 

When using AIC and SC to select models, SC tends to have larger penalty for adding terms in 
  the model.
*/





/*Chapter 12 Nonparametric model */
title 'PLot the distribution of Salary by gender';
proc univariate data=salary;
id Subj;
class Gender;
var Income;
histogram /normal;
probplot/normal(mu=est sigma=est);
run;
/*The normality assumption is not satisfied in this case, so t-test is not appropriate. 

Try Wilcoxon rank-sum test to compare income by gender. */

title 'Perform a Wilcoxon Rank-Sum Test';
proc npar1way data=Salary wilcoxon;
class Gender;
var Income;
exact;
run;
/*WILCOXON options requests wilcoxon rank sum test. If not written, other tests will be done. 
EXACT statement is not necessary, and it be used for small sample sizes only, or it will take 
  a long time to process. It can provides an exact p value except for the ones provided w/o it.
*/


/*Use univariate procedure to obtain a Wilcoxon rank test */
data difference;
set reading;
Diff=After-Before;
run;

ods select testsforlocation;
title 'Perform a Wilcoxon Signed Rank Test';
proc univariate data=difference;
var Diff;
run;
/*ods select statement restricts the output for the section called TESTS FOR LOCATION. Optional.*/


/*Perform a Kruskal-Wallis One-Way ANOVA */
/*
The Kruskal–Wallis test by ranks,is a non-parametric method for testing
whether samples originate from the same distribution.
It is used for comparing two or more independent samples of equal or different sample sizes. 
It extends the Mann–Whitney U test when there are more than two groups. 

The parametric equivalent of the Kruskal-Wallis test is the one-way analysis
of variance (ANOVA). 
A significant Kruskal-Wallis test indicates that at least one 
sample stochastically dominates one other sample.

Dunn's test,or the more powerful but less well known Conover-Iman test 
would help analyze the specific sample pairs for stochastic dominance in post hoc tests.
*/

title 'Perform a Kruskal-Wallis ANOVA';
proc npar1way data=store wilcoxon;
class Region;
var Music_Sales;
run;
/*It's the same statements with the rank test but class variables are more than two levels */

/*Compare Spread: The Ansari_Bradley Test */
/*It's an alternative to the folded F-test for comparing variances */
data twogroups;
do Group='One', 'Two';
   do Subj=1 to 15;
      input Score @;
	  output;
	end;
end;
datalines;
1 3 5 7 11 20 25 30 40 55 66 77 88 90 100
2 4 8 20 24 33 40 45 55 59 60 68 69 70 71
;
title 'Perform Ansari-Bradley Test for Spread';
proc npar1way data=twogroups ab;
class Group;
var Score;
run;
/*the single @ in the data step is an instruction to hold the line

OUTPUT statement causes an obs to be written to the resulting data set 

AB option requests Ansari-Bradley test

EXACT statement can be used if sample size is small 

AB test can be used hwen there are two groups or more, and to test for scale difference
*/


/*Convert data values into ranks */
data one;
input Subj x y;
datalines;
1 3 100
2 1 200
3 5 300
4 77 400
;
proc print data=one;
run;

proc rank data=one out=two;
var x;
ranks Rank_x;
run;

title 'Listing of data set TWO';
proc print data=two noobs;
run;
/*RANKS statment will rank the variable in VAR statement and then create a new list called rank_x

If RANKS statement is not included, the variables in VAR will be replaced by their ranks.  */

/*Run a t test on the ranks*/
title 'Run a T test';
proc rank data=salary out=rank_salary;
var Income;
ranks rank_of_salary;
run;

proc print data=rank_salary;
run;

proc ttest data=rank_salary;
class Gender;
var rank_of_salary;
run;


/*One usage of proc RANK: group the data values */
title 'Create Groups by proc RANK';
proc rank data=salary out=new_salary groups=4;
var Income;
ranks Salary_Group;
run;
proc print data=new_salary(obs=10) noobs;
run;
/*GROUPS =N divides the data values into n groups from 0 to n-1*/


/*It's often instructive to perform both a parametric and nonparametric test on the data.
 and see if they come to the similar conclusions. */







/*Chapter 13 Power and Sample Size */

/*Proc POWER or GLMPOWER */
/*Compare the sample size for an unpaired t-test */
title 'Sample Size Requirements for a T-test';
proc power;
twosamplemeans groupmeans=(20 30)(22 28) stddev=10 15 power=0.80 0.90 npergroup=.;
plot x=power min=0.7 max=0.9;
run;
/* The procedure can compute either power or sample size as long as specify by entering a 
      missing value for it*/

title 'Compute the power of an unpaired t-test';
proc power;
twosamplemeans groupmeans=20 |30 35 stddev=10 15 power=. npergroup=30 35;
plot x=n min=20 max=50;
run;


/*Compute sample size for an anova design */
title 'Compute the power for an ANOVA model';
proc power;
onewayanova groupmeans=20 |25 |30 stddev=8 10 power=.80 .90 npergroup=.;
plot x=power min=0.7 max=.90;
run;

/*Compute sample size (or power) for a difference in two proportions */
title 'Sample size for a difference in two proportions';
proc power;
twosamplefreq test=pchi groupproportions=.15 |.20 .225 .25 power=.80 .90 npergroup=.;
plot x=power min=.70 max=.90;
run;
/*Compute how many subjects needed in each group to compare a base-line proportion of 0.15 
  against three possible values: .2,.225,.25

  TEST=PCHI is Pearson's chi square, other choices include Fisher's exact test(TEST=FISHER)
    or likelihod ratio chi square test (TEST=LRCHI)*/




/*Select Random Samples */
/*SURVERYSELECT lets you take a simple random sample from a SAS data set by specifying the 
   sample size or the proportion of the input data you want to select. */

title 'Take a simple random sample';
proc surveyselect data=risk out=risk_sample method=srs samprate=0.1 seed=12434534;
run;
proc print data=risk_sample;
run;
/*SRS simple random sample method, can also be URS
  SAMPRATE: specify the rate (sample size), can also use SAMPSIZE= to specify the size directly
  SEED: 7-9 digits large number preferred*/

title 'Take a random sample with replacement';
proc surveyselect data=reading out=read_replace outhits method=urs sampsize=5 seed=12323433;
run;
proc print data=read_replace;
run;
/*OUTHITS is included so the obs which are selected more than once can be selected correctly, or
   it will produce only one single obs in the output

  URS unrestrictly random sampling */

/*Create replicate samples */
title 'Request replicate samples';
proc surveyselect data=risk out=riskrep method=srs sampsize=5 reps=3 seed=12343533;
run;
proc print data=riskrep;
run;
/*can create replicates of samples, and it's useful for bootstrap moethods and some Monte_Carlo
 techniques */


