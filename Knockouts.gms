OPTION decimals = 8
       sysout = off
       solprint = on
       reslim = 100000
       iterlim = 10000000
       domlim = 1000
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 0.0
       work = 10000000
       nlp = minos5
       mip = cplex;

$onempty

Sets
i	metabolites
/
$include AllRxns_metabs.txt
/

j	reactions
/
$include AllRxns_rxns.txt
'flavodoxin' 
/

nutrient(j)
/
$include Nutrients_M9_n2fix.txt
/

knockouts(j)
/
$include knockouts
/
;


Parameters
S(i,j)	stoichiometric matrix
/
$include AllRxns_sij.txt
*** 1.0 flxr_c -> 1.0 flxso_c + 2.0 h_c ***
'flxr[c]'.'flavodoxin' -1
'flxso[c]'.'flavodoxin' 1
'h[c]'.'flavodoxin' 1
/

rxntype(j)	reaction types
/
$include AllRxns_rxntype.txt
'flavodoxin' 1 
/

UB(j)	Upper bound on fluxes

LB(j)	Lower bound on fluxes

c(j)

maximum_biomass

flavodoxin_val

atp_val
;


*** rxntype = 1 ***
UB(j)$(rxntype(j) = 1) = 1000;
LB(j)$(rxntype(j) = 1) = 0;

*** rxntype = 0 ***
UB(j)$(rxntype(j) = 0) = 1000;
LB(j)$(rxntype(j) = 0) = -1000;

*** rxntype = 3 ***
UB(j)$(rxntype(j) = 3) = 1000;
LB(j)$(rxntype(j) = 3) = 0;

*** nutrients ***
UB(j)$(nutrient(j)) = 1000;
LB(j)$(nutrient(j)) = -1000;
LB('EX_glc(e)') = -20;
LB('EX_cbl1(e)') = -0.01;
LB('EX_o2(e)') = 0;
LB('EX_nac(e)') = -1;
LB('EX_leu_DASH_L(e)') = -1;


*** nutrients ***
*UB(j)$(inactive(j)) = 0;
*LB(j)$(inactive(j)) = 0;


***SETTING KNOCKOUTS TO 0 *****
UB(j)$(knockouts(j)) = 0;
LB(j)$(knockouts(j)) = 0;


VARIABLES
v(j)	flux of reaction j
Z	objective
;

v.lo(j)=LB(j);
v.up(j)=UB(j);
v.fx('Ec_biomass_iJO1366_core_53p95M')=0;



EQUATIONS
        obj
        massbalance(i)
        ATPmreq
        massbio
        
	;

obj..z=e=sum(j,c(j)*v(j));
massbalance(i)..sum(j,S(i,j)*v(j))=e=0;
ATPmreq..v('ATPm')=g=3.15;
massbio..v('Ec_biomass_iJO1366_WT_53p95M')=g=maximum_biomass;


c('Ec_biomass_iJO1366_WT_53p95M')=1;

***  Assemble the Model  ***
Model Biomass
/
obj
massbalance
ATPmreq
/;
Biomass.optfile=1;

Model Products
/
obj
massbalance
ATPmreq
massbio
/;
Products.optfile=1;

***  Max that Biomass eqn  ***
solve Biomass using lp maximizing z;

maximum_biomass=z.l;

file OUTPUT /autotroph_max_report_KO.txt/ ;

PUT OUTPUT;
PUT "Biomass: ";
PUT maximum_biomass:15:5/;

***Minimizing ATP***
 
c('Ec_biomass_iJO1366_WT_53p95M')=0;
c('ATPm')=1;
solve Products using lp minimizing z;

atp_val=z.l;

PUT OUTPUT;
PUT "Minimum ATP: ";
PUT atp_val:15:5/;

***Maximizing ATP***

solve Products using lp maximizing z;

atp_val=z.l;

PUT OUTPUT;
PUT "Maximum ATP: ";
PUT atp_val:15:5/;

***Minimizing Flavodoxin***
c('ATPm')=0;
c('flavodoxin')=1
solve Products using lp minimizing z;

flavodoxin_val=z.l;

PUT OUTPUT;
PUT "Minimum flavodoxin: ";
PUT flavodoxin_val:15:5/;

***Maximizing Flavodoxin***
solve Products using lp maximizing z;

flavodoxin_val=z.l;

PUT OUTPUT;
PUT "Maximum flavodoxin: ";
PUT flavodoxin_val:15:5/;





