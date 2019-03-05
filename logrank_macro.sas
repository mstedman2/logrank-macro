********************************************************************************;
*   Program:  logrank_macro.sas;
*   Purpose:  To perform a clustered logrank test;
*   input:    /usr/path/dset.sas7bdat;
*   output:   /usr/path/permuteout.sas7bdat;
*   definitions:  datain = input datset;
*                 group = variable for assigned treatment group;
*                 cluster = cluster id;
*                 censor = indicator for event=1 or censor=0;
*                 time = observation time until event or censor;
*
*  Copyright (C) 2008  Margaret Stedman;
*  
*    This program is free software: you can redistribute it and/or modify;
*    it under the terms of the GNU General Public License as published by;
*    the Free Software Foundation, either version 3 of the License, or;
*    any later version.;
*
*    This program is distributed in the hope that it will be useful,;
*    but WITHOUT ANY WARRANTY; *without even the implied warranty of;
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the;
*    GNU General Public License for more details:;
*    <http://www.gnu.org/licenses/>;
******************************************************************************;

******************************************************************************; 
******Enter directory for input dataset here:*********************************;
libname ll "N:\Projects\Other or mixed\Margaret\clustered logrank";

options nofmterr;
%macro clusterlr (datain, time, cluster, group, censor);

data one;
   set &datain;
   keep &cluster &group &time &censor;
run;

proc sql;
   select &group, count(distinct &cluster) as groupc label="Clusters", count(*) as groupn label="Observations"
   into :groupv1 - :groupv2, :groupc1 - :groupc2, :groupn1 - :groupn2
   from one
   group by &group;
quit;    

%let groups2= %eval(&groupn1 + 1);
%let totaln = %eval(&groupn1 + &groupn2);
%let totalc= %eval(&groupc1 + &groupc2);
   
proc sort data=one;
   by decending &time;
run;

data cumm_y;
    set one;
	by descending &time ;
    if _n_=1 then do;
       yy=0;
	   y_g1=0;
	   y_g2=0;
	end;
	if first.&time then do;
       dd=0;
	   d_g1=0;
	   d_g2=0;
	end;
	retain yy dd y_g1 y_g2 d_g1 d_g2;
    yy=yy+1;
	if &group=&groupv1 then y_g1=y_g1+1;
    if &group=&groupv2 then y_g2=y_g2+1;
	if &censor = 1 then do;
       dd=dd+1;
	   if &group=&groupv1 then d_g1=d_g1+1;
	   if &group=&groupv2 then d_g2=d_g2+1;
	end;
	wt=y_g1*y_g2/yy;
	time= &time;
	if last.&time then output;
	keep time yy dd y_g1 y_g2 d_g1 d_g2 wt;
run;

proc sql;
    create table teststat as
	select sum((dd * wt**2) / y_g1 / yy) as vi_1, 
	       sum((dd * wt**2) / y_g2 / yy) as vi_2, 
		   sum(d_g1 * wt / y_g1 ) as test_1,
		   sum(d_g2 * wt / y_g2 ) as test_2
	from cumm_y;

	create table newoned as
	select a.time, a.wt, a.y_g1, a.y_g2, a.yy, b.&group as group, b.&cluster as cluster, 
           b.&censor as censor, 1/a.y_g1 as y_g1_inv, 1/a.y_g2 as y_g2_inv, 1/a.yy as yy_inv 
	from cumm_y a, one b
    where a.time=b.&time;
quit;

proc sort data=newoned;
    by group cluster;
run;

proc iml;
use newoned var{time censor yy_inv y_g1 y_g2 y_g1_inv y_g2_inv wt};
read all into newone;
A=newone[1:&groupn1 , {2 8 6}];
B=newone[&groups2 : &totaln , {2 8 7}];
C= A//B;
D=C[,#];
use newoned var{cluster};
read all into csize;
%do r=1 %to &totalc;
e_pt1p=D[loc(csize=&r)][+];
e_pt1 = e_pt1//e_pt1p;
%end;
free A B C D;
A=newone[,{2 4}];
B=newone[,{2 5}];
C=A[,#];
D=B[,#];
E=repeat(C,1,&groupc1 );
F=repeat(D,1,&groupc2 );
G= E||F; 
forend=choose(G=0,0,1);   
free A B C D E F G;
A=newone[,1];
B=repeat(A,1,&totaln);
C=repeat(A`, &totaln ,1);
D=(C>=B);
%do s=1 %to &totalc;
tempp=D[1:&totaln,loc(csize=&s)][,+];
temp= temp||tempp;
%end;
free A B C D;	
A=newone[ ,{8 6 3}];
B=newone[ ,{8 7 3}];
C=A[,#];
D=B[,#];
E=repeat(C,1,&groupc1);
F=repeat(D,1,&groupc2);
fwt=E||F; 
free A B C D E F;
combp2=forend # temp # fwt;
sumpt2=combp2[+,];
e_pt2= sumpt2`;
e=e_pt1 - e_pt2;
group1=j(&groupc1 ,1,1);
group2=j(&groupc2 ,1,2);
group=group1//group2;
cluster1=1: &groupc1;
cluster2=1: &groupc2;
cluster=cluster1`//cluster2`;
fin=group || cluster || e_pt1 || e_pt2 || e;
create clustv_ijk from fin;
append from fin; 
quit;

data fin;
   set clustv_ijk;
   rename col1=group col2=cluster col3=e_pt1 col4=e_pt2 col5=e;
run;

ods listing;
proc sql;   
   create table ek as
   select group, sum(e**2) as vk
   from fin
   group by group;

   create table eachtreat as
   select sum(vk) as clustered_v
   from ek;

   title "final result";
   select (a.vi_1 + a.vi_2) as ind_var label="independent variance", b.clustered_v label="clustered variance",
  (a.test_2-a.test_1) as diff label="logrank statistic",(calculated diff/sqrt(clustered_v)) as test label="clustered logrank test",
   (1-probnorm(abs(calculated test)))*2 as p label="clustered logrank p-value"
   from teststat a, eachtreat b;

quit;

%mend;
%clusterlr (datain=ll.simdata, time=time, cluster=drid, group=treat, censor=censor);
quit;









