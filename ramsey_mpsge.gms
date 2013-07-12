$TITLE  Ramsey Model: MPSGE formulation by M.Lau, A.Pahlke and T.F.Rurherford

$ontext

SAM data at the benchmark period

           |Production| Household |  Cons.    Inv.  | Total
           |          |  L    K   |                 |
 ----------|----------|-----------|-----------------|------
 Goods     |          |           |  7.5      2.5   | 10.0
 Labor     |   7.0    |           |                 |  7.0
 Capital   |   3.0    |           |                 |  3.0
 Household |          |  7.0  3.0 |                 | 10.0
 Saving    |          |           |  2.5     -2.5   |  0.0
 ----------|----------|-----------|-----------------|------
 Total     |  10.0    |  7.0  3.0 | 10.0      0.0   |

 K0   = IO / (gamma+delta)  25  = 2.5/(0.02+0.08)
 RK0  = K0*(rbar+delta)      3   = 25*(0.04+0.08)
 L0   = Y0 - RK0            7   = 10-3
 beta = RK0 / Y0            0.3 = 3/10

$offtext

set
     t         time period             /2010*2020/
     tfirst(t) first period
     tlast(t)  last period;

tfirst(t) = yes$(ord(t) eq 1);
tlast(t)  = yes$(ord(t) eq card(t));

parameter
     gamma     growth rate             /0.02/
     rbar      Interest rate & time preference        /0.04/
     delta     depreciation rate       /0.08/
     K0        Initial Capital         /25.0/
     Y0        Benchmark income        /10.0/
     I0        Benchmark investment    /2.5/
     C0        Benchmark consumption   /7.5/
     L0        Benchmark labor input   /7.0/
     rk0       Benchmark return to capital
     L(t)      labor force
     qref(t)   Reference of quantity path
     pref(t)   Reference of price path
     beta(t)   Reference of weighted preference path
     taxk      tax rate on Capital income
;


qref(t) = (1+gamma)**(ORD(T)-1);
pref(t) = (1/(1+rbar))**(ORD(T)-1);

L(t)    = L0*qref(t);

taxk=0;



$ONTEXT

$MODEL:RAMSEY

$SECTORS:
        Y(T)    !       output
        I(T)    !       investment
        K(T)    !       capital stock
        C(T)    !       Consumption

$COMMODITIES:
        P(T)    !       output price
        RK(T)   !       return to capital
        PK(T)   !       capital price
        PL(T)   !       wage rate
        PU(T)   !       Utility price
        PKT     !       terminal capital

$CONSUMERS:
        HA      !       representative household agent

$AUXILIARY:
        TK      !       post-terminal capital stock

$PROD:Y(T) s:1
        O:P(T)          Q:Y0
        I:PL(T)         Q:L0
        I:RK(T)         Q:K0    P:(delta+rbar)  A:HA  T:taxk

$PROD:K(T)
        O:PK(T+1)       Q:(1-delta)
        O:PKT$TLAST(T)  Q:(1-delta)
        O:RK(T)         Q:1
        I:PK(T)         Q:1

$PROD:I(T)
        O:PK(T+1)       Q:1
        O:PKT$TLAST(T)  Q:1
        I:P(T)          Q:1

$PROD:C(T)
        O:PU(T)         Q:1
        I:P(T)          Q:1

$DEMAND:HA  s:1
        D:PU(T)         Q:(C0*qref(t))  P:PREF(T)
        E:PL(T)         Q:L(t)
        E:PK(TFIRST)    Q:K0
        E:PKT           Q:-1            R:TK


$CONSTRAINT:TK
        SUM(T$TLAST(T+1),  I(T+1)/I(T) - Y(T+1)/Y(T)) =E= 0;

$REPORT:
       V:QY(T)          O:P(T)       PROD:Y(t)

$OFFTEXT

$SYSINCLUDE mpsgeset RAMSEY

Y.L(T)  = Y0*qref(t);
I.L(T)  = I0*qref(t);
K.L(T)  = K0*qref(t);
C.L(T)  = C0*qref(t);

p.l(t)  = pref(t);
rk.l(t) = pref(t) * (delta+rbar);
pk.l(t) = (1+rbar) * pref(t);
pl.l(t) = pref(t);
pu.l(t) = pref(t);
pkt.l   = sum(tlast, pref(tlast));
tk.l    = k0 * (1+gamma)**card(t);

ramsey.iterlim = 0;
$include ramsey.gen
solve ramsey using mcp;

ramsey.iterlim = 1000;
$include ramsey.gen
solve ramsey using mcp;

parameter
   report(*,*,*)        report values;

set vlist /pk,Y,C,I,K/;


report(t,"pk","base") = pk.l(t)/p.l(t);
report(t,"Y","base")  = QY.l(t);
report(t,"C","base")  = C.l(t);
report(t,"I","base")  = I.l(t);
report(t,"K","base")  = K.l(t);

*------------------------------------------------
* Alternative simulation : tax on capital income
*------------------------------------------------

taxk=-0.1;

$include ramsey.gen
solve ramsey using mcp;

report(t,"pk","alt") = pk.l(t)/p.l(t);
report(t,"Y","alt")  = QY.l(t);
report(t,"C","alt")  = C.l(t);
report(t,"I","alt")  = I.l(t);
report(t,"K","alt")  = K.l(t);

report(t,vlist,"chg")=(report(t,vlist,"alt")/report(t,vlist,"base")-1)*100;

$setglobal  result   ramsey

execute_unload '%result%.gdx', report;

execute 'gdxxrw.exe %result%.gdx  par=report   rng=macro!a1    Rdim=1 Cdim=2';

execute 'del %result%.gdx';
