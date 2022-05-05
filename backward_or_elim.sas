%macro backward_or_elim(IV, DV, covariates, threshold, dataset);

/* Initialize variables */
%let cov_list=%sysfunc(compress(&covariates,(,),));
%let iteration=1;
%let num=%sysfunc(countw(&cov_list));
%let minimum=&threshold;

ods exclude all;
%do %while(&minimum<=&threshold or %eval(iteration<&num));

/*** Step 1: run full model ***/
ods output OddsRatios=OddsData_Full;
proc logistic data=&dataset descending;
class &DV / param=ref ;
model &DV = &IV &cov_list;
run;
ods output close;
proc sql noprint;
select OddsRatioEst into :fullor separated by ' '
from OddsData_Full
where Effect="%UPCASE(&IV)";
quit;

/*** Step 2: run reduced models ***/
%do i = 1 %to %sysfunc(countw(&cov_list));
%let cov_list_reduced =
%sysfunc(tranwrd(&cov_list,%scan(&cov_list,&i),));
ods output OddsRatios=OddsData_Reduced;
proc logistic data=&dataset descending;
class &DV / param=ref ;
model &DV = &IV &cov_list_reduced;
run;
data OddsData_Reduced;
length Effect $25 Removed $25;
set OddsData_Reduced;
removed="%sysfunc(scan(&cov_list,&i))";
if Effect ^= "%UPCASE(&IV)" then delete;
run;
proc append data=OddsData_Reduced base=OddsData_Merged; run;
%end;

/*** Step 3: compute effect of deleting one variable on OR ***/
data OddsData_Merged;
set OddsData_Merged;
delta = abs((OddsRatioEst-&fullor)/&fullor);
iteration=&iteration;
oddsratio=&fullor;
run;
proc sql noprint;
select min(delta) as minimum, removed into :minimum, :removedvar
from OddsData_Merged
having delta=minimum;
quit;
data OddsData_Merged;
set OddsData_Merged;
elim = "&removedvar";
run;
proc append data=OddsData_Merged base=OddsData_Final; run;
proc datasets nolist;
delete OddsData_Merged;
quit;
/*** remove &removedvar, the variable with lowest effect on OR */
%if %eval(&minimum<&threshold) %then %do;
%let cov_list=%sysfunc(tranwrd(&cov_list,&removedvar,));
%let iteration=%eval(&iteration+1);
%end;

%end;
ods exclude none;

data final;
set OddsData_Final;
run;

proc datasets nolist;
delete OddsData_Full OddsData_Reduced OddsData_Merged OddsData_Final;
quit;

ods rtf;
title "Backwards elimination procedure for independent variable &IV,
dependent variable &DV, and additional covariates at threshold &threshold";
PROC REPORT DATA=Final NOWD;
COLUMNS iteration oddsratio removed oddsratioest delta elim;
DEFINE iteration / GROUP 'Iteration';
DEFINE oddsratio / GROUP 'Full Model aOR';
DEFINE removed / GROUP 'Reduced Variable';
DEFINE oddsratioest / 'Reduced Model aOR';
DEFINE delta / FORMAT=Percent8.3 GROUP 'Change in aOR';
DEFINE elim / GROUP noprint;
break after iteration/;
compute after iteration;
line '';
endcomp;
COMPUTE delta;
IF (delta<&threshold) THEN DO;
CALL DEFINE(_col_,"STYLE","STYLE=[FONT_WEIGHT=BOLD]");
END;
ENDCOMP;
COMPUTE elim;
IF (elim = removed and delta<&threshold) THEN DO;
CALL DEFINE(_row_,"STYLE","STYLE=[BACKGROUND= cxDDDDDD]");
END;
ENDCOMP;
RUN;
ods rtf close;

%mend backward_or_elim;