
/***************************************************
*                                                  *
*               BS853 FINAL PROJECT                * 
*                                                  *
*                Log-Linear Model                  *
*                                                  *
*   Outcome - cell counts - Poisson Distribution   *
*          Variance Function - Identity            *
*              Link Funciton - Log                 *
*                                                  *
***************************************************/

dm 'odsresults; clear;';
title;

/* Data */

libname glmproj 'C:\Users\lerus\OneDrive\Documents\BU-AB\BS853-GLMs\GLM_FinalProject';

* Select for 3rd exam data (bc only exam 3 has cholesterol data);
data framingham;
set glmproj.frmgham;
if period = 3;
run;

* Creating heart disease variable by data;
data framingham2;
set framingham;
if anychd=1 or cvd=1 or PREVAP=1 or PREVMI=1 or PREVCHD=1 or PREVSTRK=1 then heart_disease=1;
else heart_disease=0;
run;
proc print data=framingham2 (obs=5);
run;

proc freq data=framingham2;
tables heart_disease;
run;

* Coding cholesterol high cholesterol and LDL variables;
data fram;
set framingham2;
if LDLC>=160 then highLDLC=1; else if LDLC<160 then highLDLC=0; else highLDLC=""; /*160 = high LDL cholesterol according to: https://medlineplus.gov/ldlthebadcholesterol.html*/
if TOTCHOL>=240 then highChol=1; else if TOTCHOL<240 then highChol=0; else highChol=""; /*>=240 cholesterol is high*/
run;
proc print data=fram (obs=10);
run;

proc freq data=fram;
tables heart_disease;
run;


/***** LDL Workflow *****/

/* Contingency Table */ 
title 'Contingency Table (LDL)';
proc freq data=fram;
table cursmoke*highLDLC*heart_disease / out=framtable1;
run;

* Odds ratio of high LDL on heart disease = 1.19;
proc freq data=fram;
table cursmoke*highLDLC*heart_disease / cmh chisq oddsratio;
run;
* Odds ratio of  heart_disease by CURSMOKE = 1;
proc freq data=fram;
table highLDLC*heart_disease*cursmoke / cmh chisq oddsratio;
run;
* Odds ratio of CURSMOKE by highLDLC = 1;
proc freq data=fram;
table heart_disease*cursmoke*highLDLC / cmh chisq oddsratio;
run;

/* MODELS */

title 'Saturated Model';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|highLDLC|heart_disease/ dist=poisson link=log obstats type3;
run;

title 'All 2-way interactions';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|highLDLC cursmoke|heart_disease highLDLC|heart_disease/ dist=poisson obstats type3;
run;


data LR1;
LR=2*(16803.1095-16801.8968);
pvalue=1-probchi(LR,1);
run;
proc print data=LR1;
title 'LRT: all 2-way vs saturated';
run;



title 'Conditional Independence of cursmoke and highLDLC';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|heart_disease highLDLC|heart_disease/ dist=poisson obstats type3;
run;


data LR2;
LR=2*(16803.1095-16801.8186);
pvalue=1-probchi(LR,2);
run;
proc print data=LR2;
title 'LRT: cursmoke|heart_disease highLDLC|heart_disease vs saturated';
run;


title 'Conditional Independence of cursmoke and heart_disease';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|highLDLC highLDLC|heart_disease/ dist=poisson obstats type3;
run;


data LR3;
LR=2*(16803.1095-16801.5480);
pvalue=1-probchi(LR,2);
run;
proc print data=LR3;
title 'LRT: cursmoke|highLDLC highLDLC|heart_disease vs saturated';
run;



title 'Conditional Independence of highLDLC and heart_disease';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|highLDLC cursmoke|heart_disease / dist=poisson obstats type3;
run;


data LR4;
LR=2*(16803.1095-16799.2979);
pvalue=1-probchi(LR,2);
run;
proc print data=LR4;
title 'LRT: cursmoke|highLDLC cursmoke|heart_disease vs saturated';
run;


title 'Joint independence of (cursmoke and highLDLC) from heart_disease';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|highLDLC heart_disease/ dist=poisson obstats type3;
run;


data LR5;
LR=2*(16803.1095-16798.9352);
pvalue=1-probchi(LR,3);
run;
proc print data=LR5;
title 'LRT: cursmoke|highLDLC heart_disease vs saturated';
run;


title 'Joint independence of (cursmoke and heart_disease) from highLDLC';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke|heart_disease highLDLC / dist=poisson obstats type3;
run;


data LR6;
LR=2*(16803.1095-16799.2058);
pvalue=1-probchi(LR,3);
run;
proc print data=LR6;
title 'LRT: cursmoke|heart_disease highLDLC vs saturated';
run;


title 'Joint independence of (highLDLC and heart_disease) from cursmoke';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = highLDLC|heart_disease cursmoke / dist=poisson obstats type3;
run;


data LR7;
LR=2*(16803.1095-16801.4559);
pvalue=1-probchi(LR,3);
run;
proc print data=LR7;
title 'LRT: highLDLC|heart_disease cursmoke vs saturated';
run;


title 'Mutual independence model';
ods select ModelFit;
proc genmod data=framtable1;
class cursmoke highLDLC heart_disease;
model count = cursmoke highLDLC heart_disease / dist=poisson obstats type3;
run;


data LR8;
LR=2*(16803.1095-16798.8431);
pvalue=1-probchi(LR,4);
run;
proc print data=LR8;
title 'LRT: mutual independence vs saturated';
run;


title 'Intercept only model';
ods select ModelFit;
proc genmod data=framtable1;
model count = / dist=poisson obstats type3;
run;


data LR9;
LR=2*(16803.1095-16350.7649); * Mutual indep LL: 16798.8431, (LDL|HD, Smoke) LL: 16801.4559;
pvalue=1-probchi(LR,1);
run;
proc print data=LR9;
title 'LRT: mutual independence vs (LDL|HD, Smoke)';
run;


/* FINAL MODEL */

title 'Joint independence of (highLDLC and heart_disease) from cursmoke';
proc genmod data=framtable1;
class cursmoke(ref='0') highLDLC(ref='0') heart_disease(ref='0');
model count = highLDLC|heart_disease cursmoke / dist=poisson obstats type3;
run;


/**********************************************************************/


/***** TotChol Workflow *****/

/* Contingency Table */ 
title 'Contingency Table (TotChol)';
proc freq data=fram;
table cursmoke*highChol*heart_disease / out=framtable2;
run;

* Odds ratio of CURSMOKE by highChol = 0.840;
proc freq data=fram;
table heart_disease*cursmoke*highChol / cmh chisq oddsratio;
run;
* Odds ratio of highChol by heart_disease = 1;
proc freq data=fram;
table cursmoke*highChol*heart_disease / cmh chisq oddsratio;
run;
* Odds ratio of heart_disease by CURSMOKE = 1;
proc freq data=fram;
table highChol*heart_disease*cursmoke / cmh chisq oddsratio;
run;


/* MODELS */

title 'Saturated Model';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|highChol|heart_disease / dist=poisson obstats type3;
run;

title 'All 2-way interactions';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|highChol cursmoke|heart_disease highChol|heart_disease/ dist=poisson obstats type3;
run;

data LR21;
LR=2*(16801.1717-16799.8689);
pvalue=1-probchi(LR,1);
run;
proc print data=LR21;
title 'LRT: all 2-way vs saturated';
run;


title 'Conditional Independence of cursmoke and highChol';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|heart_disease highChol|heart_disease/ dist=poisson obstats type3;
run;


data LR22;
LR=2*(16801.1717-16797.1639);
pvalue=1-probchi(LR,2);
run;
proc print data=LR22;
title 'LRT: cursmoke|heart_disease highChol|heart_disease vs saturated';
run;


title 'Conditional Independence of cursmoke and heart_disease';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|highChol highChol|heart_disease/ dist=poisson obstats type3;
run;


data LR23;
LR=2*(16801.1717-16799.5507);
pvalue=1-probchi(LR,2);
run;
proc print data=LR23;
title 'LRT: cursmoke|highChol highChol|heart_disease vs saturated';
run;



title 'Conditional Independence of highChol and heart_disease';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|highChol cursmoke|heart_disease / dist=poisson obstats type3;
run;

data LR24;
LR=2*(16801.1717-16799.0168);
pvalue=1-probchi(LR,2);
run;
proc print data=LR24;
title 'LRT: cursmoke|highChol cursmoke|heart_disease vs saturated';
run;

title 'Joint independence of (cursmoke and highChol) from heart_disease';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|highChol heart_disease/ dist=poisson obstats type3;
run;


data LR25;
LR=2*(16801.1717-16798.6542);
pvalue=1-probchi(LR,3);
run;
proc print data=LR25;
title 'LRT: cursmoke|highChol heart_disease vs saturated';
run;

data LR25a;
LR=2*(16799.0168-16798.6542);
pvalue=1-probchi(LR,1);
run;
proc print data=LR25a;
title 'LRT: cursmoke|highChol heart_disease vs conditional on smoke';
run;

data LR25b;
LR=2*(16799.5507-16798.6542);
pvalue=1-probchi(LR,1);
run;
proc print data=LR25b;
title 'LRT: cursmoke|highChol heart_disease vs conditional on highChol';
run;


title 'Joint independence of (cursmoke and heart_disease) from highChol';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke|heart_disease highChol / dist=poisson obstats type3;
run;


data LR26;
LR=2*(16801.1717-16796.2673);
pvalue=1-probchi(LR,3);
run;
proc print data=LR26;
title 'LRT: cursmoke|heart_disease highChol vs saturated';
run;


title 'Joint independence of (highChol and heart_disease) from cursmoke';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = highChol|heart_disease cursmoke / dist=poisson obstats type3;
run;

data LR27;
LR=2*(16801.1717-16796.8012);
pvalue=1-probchi(LR,3);
run;
proc print data=LR27;
title 'LRT: highChol|heart_disease cursmoke vs saturated';
run;


title 'Mutual independence model';
ods select ModelFit;
proc genmod data=framtable2;
class cursmoke highChol heart_disease;
model count = cursmoke highChol heart_disease / dist=poisson obstats type3;
run;


data LR28;
LR=2*(16801.1717-16795.9047);
pvalue=1-probchi(LR,4);
run;
proc print data=LR28;
title 'LRT: mutual independence vs saturated';
run;


title 'Intercept only model';
ods select ModelFit;
proc genmod data=framtable2;
model count = / dist=poisson obstats type3;
run;

data LR29;
LR=2*(16801.1717-16350.7649);
pvalue=1-probchi(LR,7);
run;
proc print data=LR29;
title 'LRT: intercept only vs saturated';
run;


/* FINAL MODEL */

title 'Joint independence of (cursmoke and highChol) from heart_disease';
proc genmod data=framtable2;
class cursmoke(ref='0') highChol(ref='0') heart_disease(ref='0');
model count = cursmoke|highChol heart_disease / dist=poisson obstats type3;
run;


