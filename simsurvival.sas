*********************************************************************************;
*   program: simsurvival.sas;
*   author:  margaret stedman;
*   date:    Jan 19, 2009;
*   purpose: to simulate survival data to demonstrate logrank macro ;
********************************************************************************;

libname ll "N:\Projects\Other or mixed\Margaret\clustered logrank";

%macro sim_iml(r00,r01,gamma,ndrs,ntreat,npats,seed,totpats,nsim,name);

data _null_;
   shape=1 / &gamma;
   b00=(-1) * &r00 * &gamma;
   b01=(-1) * &r01 * &gamma;
   call symput('shape',shape);
   call symput('b00',b00);
   call symput('b01',b01);
run;

proc iml;
dref=j(&nsim,&ndrs,.); 
time=j(&nsim,&totpats,.);
cens=j(&nsim,&totpats,1);
tout=j(1,1,.);
lparm=j(1,1,.);
truetreat=j(1,&totpats,&r01); 
shape=j(1,&totpats, &gamma);
treatgrp=repeat(0:1,1,&ntreat);
drid=repeat(1:&ndrs,1,&npats); 
call randseed (&seed); 
call randgen(dref,'NORMAL'); 
mod=repeat(dref,1,&npats) + repeat(treatgrp,&nsim,&npats) # &b01 + &b00; 
do f=1 to &nsim;
	do g=1 to &totpats;
	   lparm=mod[f,g];
	   lambda=exp(lparm);
	   call randgen(tout,'WEIB', &shape ,(lambda**(-1 * &shape )));
	   time[f,g]=tout;
	   if tout>296 then do;
         cens[f,g]=0;
         time[f,g]=296;
	   end;
	end;
end;
treat=repeat(treatgrp,1,&npats); 
fin=truetreat`||shape`||drid`||treat`||time`||cens`;
create fin&name from fin;
append from fin;
quit;
%mend sim_iml;
%sim_iml(r00=7.6434, r01=-.10, gamma=1.1, ndrs=458, ntreat=229, npats=8, seed=345, totpats=3664, nsim=1, name=1);

data ll.simdata;
  set fin1;
  rename col1=truetreat col2=shape col3=drid col4=treat col5=time col6=censor;
run;
