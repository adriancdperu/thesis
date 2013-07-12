$title  Static General Equilibrium for Japan based on 2005 IO table

option solprint=off;

*====================================================================
*   Set definitions for 4 x 4 model
*====================================================================

SET
  i      Sectors
             /AGR    Agriculture
              MFG    Manufacturing
              SRV    Services
              ENE    Energy/;
SET
  en(i)  Energy Sector  /ene/;


alias (i,j);

*====================================================================
*   SAM data of the Japanese economy based on 2005 IO table
*====================================================================

* set file name of your IO data

$setglobal    iodata     io4x4

$setglobal    result   %iodata%_dyn_co2

parameter  IO_data(*,*),CO2_data(*,*);

$call gdxxrw   %iodata%.xlsx  par=IO_data   rng=io!A1:L9  Rdim=1 Cdim=1
$gdxin %iodata%.gdx
$load  IO_data
$gdxin

$call gdxxrw   %iodata%.xlsx  par=CO2_data  rng=co2!A1:F2  Rdim=1 Cdim=1
$gdxin %iodata%.gdx
$load  CO2_data
$gdxin
execute 'del %iodata%.gdx';

display io_data, co2_data;

*====================================================================
*   Variables definitions
*====================================================================

PARAMETERS
  Y0(i)        Domestic output by commodities
  YI0(i)       Domestic supply for domestic use by commodities
  YA0(i)       Armington supply by commodities
  IJ0(i,j)     Intermediate inputs
  CI0(i)       Private consumption by commodities
  C0           Total private consumption
  GI0(i)       Government consumption by commodities
  G0           Total government consumption
  INVI0(i)     Investment by commodities
  INV0         Total investment
  XI0(i)       Export by commodities
  X0           Total export
  MI0(i)       Import by commodities
  M0           Total import
  LI0(i)       Labor input by industry
  L0           Total labor input
  KI0(i)       Capial input by industry
  K0           Total capital input
  LES0         Leisure
  UTIL0        Utility
*
  HSAVE0       Household saving
  GSAVE0       Government saving
  BOPSURP0     Balance of payment surplus
*
  TAXM0(i)     Import tariff of commodity i
*
  txy0(i)      Tax rate on output of sector i
  txm0(i)      Tax rate on import of commodity i
  txl0         Tax rate on labor income
  txk0         Tax rate on capital income
*
  pl0          Wage
  rk0          Rate of return to capital
  py0(i)       Price of output  of sector i
  pm0(i)       Import price (tariff included) of commodity i
*
;

*  Intermediate inputs

IJ0(i,j)  = IO_data(i,j);

*  Private consumption by commodities

CI0(i)    = io_data(i,"cp");
C0        = sum(i,CI0(i));

*  Government consumption by commodities

GI0(i)    = io_data(i,"gov");
G0        = sum(i,GI0(i));

*  Investment by commodities

INVI0(i)  = io_data(i,"inv");
INV0      = sum(i,INVI0(i));

*  Export by commodities

XI0(i)    = io_data(i,"export");
X0        = sum(i,XI0(i));

*  Import by commodities

MI0(i)    = -io_data(i,"import");
M0        = sum(i,MI0(i));
TAXM0(i)  = -io_data(i,"tariff");

*  Domestic output

Y0(i)     = sum(j,IJ0(j,i)) + io_data("cap",i) + io_data("lab",i)  + io_data("taxi",i);

*  Domestic supply for domestic output

YI0(i)    = Y0(i) -XI0(i);

*  Armington demand

YA0(i)    = YI0(i)+MI0(i)+TAXM0(i);

* Basic data for SAM matrix

Parameter
   VTAXK0   Capital Income tax
   VTAXL0   Labor income tax;

VTAXK0 = 16.771;
VTAXL0 = 24.578;


*  Capital income and tax

KI0(j)    = IO_data("cap",j);
txk0      = VTAXK0/(sum(j,KI0(j))-VTAXK0);
KI0(i)    = io_data("cap",i)/(1+txk0);
K0        = sum(i,KI0(i));
rk0       = 1 + txk0;

*  Labor income and tax

LI0(j)    = io_data("lab",j);
txl0      = VTAXL0/(sum(j,LI0(j))-VTAXL0);
LI0(i)    = io_data("lab",i)/(1+txl0);
L0        = sum(i,LI0(i));
pl0       = 1 + txl0;


* Ouput tax

txy0(i)   = io_data("taxi",i)/Y0(i);
py0(i)    = 1 - txy0(i);

* Import tax

txm0(i)$io_data(i,"import")   = io_data(i,"tariff")/io_data(i,"import");
pm0(i)    = 1 + txm0(i);

display txy0,txm0,txk0,txl0;

LES0      =   75;
UTIL0     =   C0 + LES0;
HSAVE0    =   K0 + L0 - C0;
GSAVE0    =   VTAXK0 + VTAXL0
            + sum(j,IO_data("taxi",j))
            - sum(i,IO_data(i,"tariff")) - G0;
BOPSURP0  =   X0 - M0;

Parameter  CHK;

CHK = HSAVE0 + GSAVE0 - BOPSURP0 - INV0;

display LES0, UTIL0, HSAVE0, GSAVE0, BOPSURP0, CHK;

*===================================================================
*  Dynamic Parameter
*===================================================================

SCALARS
  gamma        Population and productivity growth rate          /0.02/,
  delta        Depreciation rate                                /0.08/,
  rbar         Interest rate                                    /0.04/,
  esubt        Intertemporal elasticity of substitution         /0.5/;


SETS
  t            period  /2005*2020/,
  tfirst(t)    First period in the model,
  tlast(t)     Last period in the model;

tfirst(t) = yes$(ord(t) eq 1);
tlast(t)  = yes$(ord(t) eq card(t));


PARAMETERS
  pref(t)      Reference price path (present value price index),
  qref(t)      Reference quantity path (index of population size);

qref(t)    = (1 + gamma)**(ord(t)-1);
pref(t)    = 1 / (1 + rbar)**(ord(t)-1);

*===================================================================
*      Calibration of initial steady state capital stock
*===================================================================

PARAMETERS
  PR0          Base year user cost of capital
  PK0          Base year price of capital
  R0           Base year capital income (Net of tax)
  CADEF0       Base year current account deficit
  GDEF0        Base year government deficit
  LSINC0       Base year total household lumpsum receipts
  CI0_0(i)     Consumption not adjusted
  C0_0         Sum of consumption not adjusted
  INVI0_0(i)   Investment not adjusted
  INV0_0       Sum of investment not adjusted
  KT0          Terminal capital
;

PR0      = (rbar + delta);
PK0      = (1 + rbar);
R0       = K0;

LSINC0   = 0;

CADEF0   = bopsurp0;
GDEF0    = gsave0;

*  Initial steady state capital stock

K0       = R0/(rbar + delta);

CI0_0(i)   = CI0(i);
C0_0       = C0;
INVI0_0(i) = INVI0(i);
INV0_0     = INV0;

*  Adjusting investment and consumption to be constistent with steady state

INV0       = (delta + gamma)*K0;
INVI0(i)   = INVI0(i)*INV0/INV0_0;
C0 = C0_0 + (INV0_0 - INV0);
CI0(i) = CI0_0(i) + ( INVI0_0(i) - INVI0(i));
display hsave0, C0_0, C0, INV0_0, INV0, CI0_0, CI0, INVI0_0, INVI0;

UTIL0     = C0+LES0;
HSAVE0    = hsave0 - (INV0_0 - INV0);

display util0,hsave0


*===================================================================
*   Tax parameters
*===================================================================

parameters
   txy(i,t)    Tax rate on output
   txl(t)      Tax rate on labor income
   txk(t)      Tax rate on capital income
   txm(i,t)    Tax rate on import
;

txy(i,t)   = txy0(i);
txl(t)     = txl0;
txk(t)     = txk0;
txm(i,t)   = txm0(i);


parameter
   co2lim(t)   limit of co2 emission;

co2lim(t)=sum(j,CO2_data("ene",j))*2;


*===================================================================
*   Dynamic MPSGE model
*===================================================================

$ONTEXT

$MODEL:Japan

$SECTORS:
  Y(i,t)           ! Domestic production
  YA(i,t)          ! Armington aggregation of imports and domestic use
  MI(i,t)$MI0(i)   ! Import
  XI(i,t)$XI0(i)   ! Export
  GOV(t)           ! Government demand
  INV(t)           ! Gross investment demand
  KS(t)            ! Capital stock
  LS(t)            ! Labor supply
  UTIL(t)          ! Private utility

$COMMODITIES:
  PA(i,t)          ! Price of Armington aggregates
  PD(i,t)          ! Price of domestic market price
  PM(i,t)$MI0(i)   ! Price of import
  PX(i,t)$XI0(i)   ! Price of export
  PFX(t)           ! Price of foreign exchange
  PGOV(t)          ! Price of government expendtitures
  PLS(t)           ! Price of Labor
  PL(t)            ! Price of Leisure
  RK(t)            ! Rate of return to capital
  PK(t)            ! Price of Capital
  PUTIL(t)         ! Price of private consumption
  PCO2(t)          ! CO2 price
  PKT              ! Post terminal capiatl constraint

$CONSUMERS:
  HA               ! Household agent
  GOVT(t)          ! Government (tax collector)

$AUXILIARY:
        TK         ! Post-terminal capital stock

*====================================================================
*   CRTS production with a CET between domestic outputs and exports:
*====================================================================

$PROD:Y(i,t)            t:2   s:0   ve(s):0.1  vc(ve):0  va(ve):1
  O:PD(i,t)             Q:YI0(i)       P:py0(i)    A:GOVT(t)   T:txy(i,t)
  O:PX(i,t)$XI0(i)      Q:XI0(i)       P:py0(i)    A:GOVT(t)   T:txy(i,t)
  I:PA(j,t)$(not en(j)) Q:IJ0(j,i)
  I:PA(j,t)$en(j)       Q:IJ0(j,i)                                  vc:
  I:PCO2(t)             Q:CO2_data("ene",i)  P:1.0e-6               vc:
  I:RK(t)               Q:KI0(i)       p:rk0       A:GOVT(t)   T:txk(t)  va:
  I:PLS(t)              Q:(LI0(i)*PL0)                              va:

*====================================================================
*   CRTS production with a CET between domestic outputs and exports:
*====================================================================

$PROD:Y(i,t)$(not XI0(i))     t:2   s:0   ve(s):0.1  vc(ve):0  va(ve):1
  O:PD(i,t)             Q:YI0(i)       P:py0(i)    A:GOVT(t)   T:txy(i,t)
  I:PA(j,t)$(not en(j)) Q:IJ0(j,i)
  I:PA(j,t)$en(j)       Q:IJ0(j,i)                                  vc:
  I:PCO2(t)             Q:CO2_data("ene",i)  P:1.0e-6               vc:
  I:RK(t)               Q:KI0(i)       p:rk0       A:GOVT(t)   T:txk(t)  va:
  I:PLS(t)              Q:(LI0(i)*PL0)                              va:


*====================================================================
*   Armington Demand
*====================================================================

$PROD:YA(i,t)    s:2
  O:PA(i,t)             Q:YA0(i)
  I:PD(i,t)             Q:YI0(i)
  I:PM(i,t)             Q:(pm0(i)*MI0(i))

*====================================================================
*   Imports
*====================================================================

$PROD:MI(i,t)$MI0(i)
  O:PM(i,t)             Q:(pm0(i)*MI0(i))
  I:PFX(t)              Q:MI0(i)        P:pm0(i)    A:GOVT(t)   T:txm(i,t)

*====================================================================
*   Exports
*====================================================================

$PROD:XI(i,t)$XI0(i)
  O:PFX(t)              Q:XI0(i)
  I:PX(i,t)             Q:XI0(i)

*====================================================================
*   Capital stock
*====================================================================

$PROD:KS(t)
  O:PKT$TLAST(t)        Q:(K0*(1-DELTA))
  O:PK(t+1)             Q:(K0*(1-DELTA))
  O:RK(t)               Q:R0
  I:PK(t)               Q:K0

*====================================================================
*   Investment:
*====================================================================

$PROD:INV(t)
  O:PKT$TLAST(t)        Q:INV0
  O:PK(t+1)             Q:INV0
  I:PA(i,t)             Q:INVI0(i)

*====================================================================
*   Labor supply
*====================================================================

$PROD:LS(t)
  O:PLS(t)              Q:(L0*PL0)
  I:PL(t)               Q:L0              P:PL0      A:GOVT(t)   T:txl(t)

*====================================================================
*   Household utility with consumption
*====================================================================

$PROD:UTIL(t)   s:1  en(s):0
  O:PUTIL(t)             Q:UTIL0
  I:PA(i,t)$(not en(i))  Q:CI0(i)
  I:PA(i,t)$en(i)        Q:CI0(i)                            en:
  I:PCO2(t)              Q:CO2_data("ene","cp")  P:1.0e-6    en:
  I:PL(t)                Q:LES0

*====================================================================
*   Government demands
*====================================================================

$PROD:GOV(t)  s:1
  O:PGOV(t)             Q:G0
  I:PA(i,t)             Q:GI0(i)

*====================================================================
*   Household agents:
*====================================================================

$DEMAND:HA  s:esubt
  D:PUTIL(t)            Q:(UTIL0*qref(t))  P:pref(t)
  E:PL(t)               Q:(UTIL0*qref(t))
  E:PK(tfirst)          Q:K0
  E:PKT                 Q:(-1)          R:TK
  E:PGOV(t)             Q:(GSAVE0*qref(t))
  E:PFX(t)              Q:(-BOPSURP0*qref(t))
  E:PCO2(t)             Q:co2lim(t)

*====================================================================
*   Government agents:
*====================================================================

$DEMAND:GOVT(t)
  D:PGOV(t)             Q:G0
  E:PGOV(t)             Q:(-GSAVE0)

*====================================================================
*   Terminal constraints:
*====================================================================

$CONSTRAINT:TK
  SUM(T$TLAST(T), INV(T)/INV(t-1) - UTIL(T)/UTIL(t-1)) =E= 0;

*====================================================================
*   Report:
*====================================================================

$REPORT:
  V:QD(i,t)            O:PD(i,t)      PROD:Y(i,t)
  V:QX(i,t)            O:PX(i,t)      PROD:Y(i,t)
  V:K(i,t)             I:RK(t)        PROD:Y(i,t)
  V:L(i,t)             I:PLS(t)       PROD:Y(i,t)
  V:LSUP(t)            I:PL(t)        PROD:LS(t)
  V:LES(t)             I:PL(t)        PROD:UTIL(t)
  V:CP(i,t)            I:PA(i,t)      PROD:UTIL(t)
  V:G(i,t)             I:PA(i,t)      PROD:GOV(t)
  V:INVEST(i,t)        I:PA(i,t)      PROD:INV(t)
  V:EX(i,t)$Xi0(i)     I:PX(i,t)      PROD:XI(i,t)
  V:IM(i,t)$MI0(i)     O:PM(i,t)      PROD:MI(i,t)
  V:VUTIL(t)           O:PUTIL(t)     PROD:UTIL(t)
  V:CO2_I(i,t)         I:PCO2(t)      PROD:Y(i,t)
  V:CO2_h(t)           I:PCO2(t)      PROD:UTIL(t)

$OFFTEXT

$SYSINCLUDE mpsgeset Japan

*====================================================================
*   Standard solution:
*====================================================================
Y.l(i,t)    = qref(t);
YA.L(i,t)   = qref(t);
MI.L(i,t)   = qref(t);
XI.L(i,t)   = qref(t);
KS.L(t)     = qref(t);
LS.l(t)     = qref(t);
INV.L(t)    = qref(t);
UTIL.L(t)   = qref(t);
GOV.L(t)    = qref(t);
TK.L        = K0*(1+gamma)**CARD(t);

PA.L(i,t)   = pref(t);
PD.L(i,t)   = pref(t);
PM.L(i,t)   = pref(t);
PX.L(i,t)   = pref(t);
PFX.L(t)    = pref(t);
PK.L(t)     = (1+rbar)*pref(t);
PKT.L       = SUM(TLAST,PK.L(TLAST)/(1+rbar));
RK.L(t)     = pref(t);
PL.L(t)     = pref(t);
PLS.l(t)    = pref(t);
PUTIL.L(t)  = pref(t);
PGOV.L(t)   = pref(t);

PCO2.L(t)   = 0;

japan.iterlim=0;
$include Japan.gen
solve Japan using mcp;

japan.iterlim=1000;
$include Japan.gen
solve Japan using mcp;

parameters report(*,*,*), report_i(*,*,*,*);

report_i(t,"Y",i,"base")      = QD.L(i,t)+QX.L(i,t);
report_i(t,"K",i,"base")      = K.L(i,t);
report_i(t,"L",i,"base")      = L.L(i,t);
report_i(t,"PD",i,"base")     = PD.L(i,t);
report_i(t,"PX",i,"base")     = PX.L(i,t);
report_i(t,"PM",i,"base")     = PM.L(i,t);
report_i(t,"PA",i,"base")     = PA.L(i,t);
report_i(t,"CP",i,"base")     = CP.L(i,t);
report_i(t,"G",i,"base")      = G.L(i,t);
report_i(t,"INV",i,"base")    = INVEST.L(i,t);
report_i(t,"EX",i,"base")     = EX.L(i,t);
report_i(t,"IM",i,"base")     = IM.L(i,t);

report(t,"CP","base")     = sum(i,CP.L(i,t));
report(t,"G","base")      = sum(i,G.L(i,t));
report(t,"INV","base")    = sum(i,INVEST.L(i,t));
report(t,"EX","base")     = sum(i,EX.L(i,t));
report(t,"IM","base")     = sum(i,IM.L(i,t));
report(t,"GDP","base")    = report(t,"CP","base") + report(t,"G","base")
                        + report(t,"INV","base")
                        + report(t,"EX","base") - report(t,"IM","base");
report(t,"LS","base")     = LSUP.L(t);
report(t,"LES","base")    = LES.L(t);
report(t,"RK","base")     = RK.L(t);
report(t,"PLS","base")    = PLS.L(t);
report(t,"PL","base")     = PL.L(t);
report(t,"PUTIL","base" ) = PUTIL.L(t);
report(t,"UTIL","base")   = UTIL.L(t);
report(t,"VUTIL","base")  = VUTIL.L(t);
report(t,"CO2","base")    = sum(i,CO2_I.L(i,t))+CO2_h.L(t);
report(t,"PCO2","base")   = PCO2.L(t)*1000000;


*====================================================================
*   Counter factural Solution
*====================================================================

* Import tariff zero

*txm(i,t)=0;

* Increase in output tax rate

*txy(i,t) = txy(i,t)+0.10;

* Increase in labor income

*txl(t) = txl(t) + 0.10;

* Increase in capital income

*txk(t) = txk(t) - 0.10;

* Increase in energy tax

*txm("ene",t)=txm("ene",t)+0.10;
*txy("ene",t)=txy("ene",t)+0.04;


CO2lim("2005")= 1181;
CO2lim("2006")= 1153;
CO2lim("2007")= 1125;
CO2lim("2008")= 1098;
CO2lim("2009")= 1072;
CO2lim("2010")= 1046;
CO2lim("2011")= 1022;
CO2lim("2012")=  997;
CO2lim("2013")=  973;
CO2lim("2014")=  950;
CO2lim("2015")=  927;
CO2lim("2016")=  905;
CO2lim("2017")=  884;
CO2lim("2018")=  863;
CO2lim("2019")=  842;
CO2lim("2020")=  822;


$include Japan.gen
solve Japan using mcp;

report_i(t,"Y",i,"alt")      = QD.L(i,t)+QX.L(i,t);
report_i(t,"K",i,"alt")      = K.L(i,t);
report_i(t,"L",i,"alt")      = L.L(i,t);
report_i(t,"PD",i,"alt")     = PD.L(i,t);
report_i(t,"PX",i,"alt")     = PX.L(i,t);
report_i(t,"PM",i,"alt")     = PM.L(i,t);
report_i(t,"PA",i,"alt")     = PA.L(i,t);
report_i(t,"CP",i,"alt")     = CP.L(i,t);
report_i(t,"G",i,"alt")      = G.L(i,t);
report_i(t,"INV",i,"alt")    = INVEST.L(i,t);
report_i(t,"EX",i,"alt")     = EX.L(i,t);
report_i(t,"IM",i,"alt")     = IM.L(i,t);

report(t,"CP","alt")     = sum(i,CP.L(i,t));
report(t,"G","alt")      = sum(i,G.L(i,t));
report(t,"INV","alt")    = sum(i,INVEST.L(i,t));
report(t,"EX","alt")     = sum(i,EX.L(i,t));
report(t,"IM","alt")     = sum(i,IM.L(i,t));
report(t,"GDP","alt")    = report(t,"CP","alt") + report(t,"G","alt")
                        + report(t,"INV","alt")
                        + report(t,"EX","alt") - report(t,"IM","alt");
report(t,"LS","alt")     = LSUP.L(t);
report(t,"LES","alt")    = LES.L(t);
report(t,"RK","alt")     = RK.L(t);
report(t,"PLS","alt")    = PLS.L(t);
report(t,"PL","alt")     = PL.L(t);
report(t,"PUTIL","alt" ) = PUTIL.L(t);
report(t,"UTIL","alt")   = UTIL.L(t);
report(t,"VUTIL","alt")  = VUTIL.L(t);
report(t,"CO2","alt")    = sum(i,CO2_I.L(i,t))+CO2_h.L(t);
report(t,"PCO2","alt")   = PCO2.L(t)*1000000;


*====================================================================
*   Ratio of Change
*====================================================================

set var_i
       /Y,K,L,PD,PX,PM,PA,CP,G,INV,EX,IM/
;
set var
       /CP,G,INV,EX,IM,GDP,LS,LES,RK,PLS,PL,PUTIL,UTIL,VUTIL,CO2,PCO2/
;

report_i(t,var_i,i,"chg")$report_i(t,var_i,i,"base")
               = (report_i(t,var_i,i,"alt")/report_i(t,var_i,i,"base")-1)*100;
report(t,var,"chg")$report(t,var,"base")
               = (report(t,var,"alt")/report(t,var,"base")-1)*100;


*====================================================================
*   Display simulation results
*====================================================================

display report, report_i;


execute_unload '%result%.gdx', report_i,report;

execute 'gdxxrw.exe %result%.gdx  par=report_i rng=industry!a1  Rdim=1 Cdim=3';
execute 'gdxxrw.exe %result%.gdx  par=report   rng=macro!a1    Rdim=1 Cdim=2';

execute 'del %result%.gdx';
