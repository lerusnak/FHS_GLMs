
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

data fram;
set glmproj.frmgham;
run;

data fram1;
set fram;
where period = 3;
if LDLC>=160 then highLDLC=1; else highLDLC=0; /*160 = high LDL cholesterol according to: https://medlineplus.gov/ldlthebadcholesterol.html*/
run;
proc print data=fram1 (obs=10);
run;


/* Contingency Table */ 

proc freq data=fram1;
table cursmoke*highLDLC*anychd / out=framtable;
run;


/* MODELS */

title 'Saturated Model';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|highLDLC|anychd / dist=poisson obstats type3;
run;

title 'All 2-way interactions';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|highLDLC cursmoke|anychd highLDLC|anychd/ dist=poisson obstats type3;
run;

title 'Conditional Independence of cursmoke and highLDLC';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|anychd highLDLC|anychd/ dist=poisson obstats type3;
run;

title 'Conditional Independence of cursmoke and anychd';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|highLDLC highLDLC|anychd/ dist=poisson obstats type3;
run;

title 'Conditional Independence of highLDLC and anychd';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|highLDLC cursmoke|anychd / dist=poisson obstats type3;
run;

title 'Joint independence of (cursmoke and highLDLC) from anychd';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|highLDLC anychd/ dist=poisson obstats type3;
run;

title 'Joint independence of (cursmoke and anychd) from highLDLC';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke|anychd highLDLC / dist=poisson obstats type3;
run;

title 'Joint independence of (highLDLC and anychd) from cursmoke';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = highLDLC|anychd cursmoke / dist=poisson obstats type3;
run;

title 'Mutual independence model';
ods select ModelFit;
proc genmod data=framtable;
class cursmoke highLDLC anychd;
model count = cursmoke highLDLC anychd / dist=poisson obstats type3;
run;
