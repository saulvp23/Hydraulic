$TITLE Pump scheduling smallest
$ontext
version: 2.0
author of this version: sofdem@gmail.com
characteristics: 9-period model for free GAMS version
$offtext

*$offlisting
* stops the echo print of the input file
*$offsymxref
* stops the print of a complete cros-reference list of symbols
*option limcol = 0;
* stops the print of the column listing
option limrow = 0;
* stops the print of the equation listing
*option solprint = off;
* stop report solution

Sets
     n          nodes      / s, j1, j2, r1, r2, r3, r4 /
     j(n)       junctions  / j1, j2 /
     r(n)       reservoirs / r1, r2, r3, r4 /
     l(n,n)     pipes      / s.j1, j1.j2, j1.r1, j1.r4, j2.r2, j2.r3 /
     t          1-hour periods / t1*t9 /
     night(t)   night periods  / t1*t4 /
     c          pump class   / small /
     d          pump number  / p1*p3 /
     k(c,d)     pump type    / small.p1*p3 /

* compatibility with gams < 24.2
*     kp(c,d)    pump type bis / small.p1*p3 /
     degree     polynomial degrees / 0*2 /
     tfirst(t);

     alias (n,np)
     alias (k,kp)
     alias (d,dp);
Scalar
     height0      reference height at the source (m)      / 0 /
     tariffnight  electricity hourly tariff at night (euro.kWh^-1)  / 0.02916 /
     tariffday    electricity hourly tariff at day (euro.kWh^-1)    / 0.04609 /;

Parameter tariff(t)   electricity tariff;
    tariff(t)        = tariffday;
    tariff(night(t)) = tariffnight;

Parameters

    height(n)   elevation at each node relative to s (m)
                / s 0, j1 30, j2 30, r1 50, r2 50, r3 45, r4 35 /

    surface(r)  mean surface of each reservoir (m^2)
                / r1 80, r2 80, r3 80, r4 80 /

    vmin(r)     minimal volume of each reservoir (m^3)
                / r1 100, r2 100, r3 100, r4 100 /

    vmax(r)     maximal volume of each reservoir (m^3)
                / r1 300, r2 300, r3 300, r4 300 /  ;

Parameter vinit(r) initial volume of each reservoir ;
    vinit(r) = vmin(r);

tfirst(t)  =  yes$(ord(t) = 1);

Parameter
         Qmax(c) Qmax en (m^3.h^-1)
                 / small 100 /
* a polynomial is represented as the list of coefficients for each term degree
Table psi(c,degree) quadratic fit of the service pressure (m) div( flow (m^3.h^-1) for each class of pumps)
                  0              1            2
      small       63.0796        0            -0.0064085;

Table gamma(c,degree) linear fit of the electrical power (euro.kWh^-1) div( flow (m^3.h^-1) for each class of pumps )
                  0              1
      small       3.81101     0.09627;

Table demand(r,t) demand in water at each reservoir each hour (m^3.h^-1)
     t1     t2    t3    t4     t5     t6    t7      t8      t9
r1   9.83   5.0   3.67  6.5    5.67   7.5   3.0     3.0     2.0
r2   44.83  18.0  0.0   0.0    0.0    0.0   0.0     45.0    51.67
r3   14.0   13.33 25.5  11.0   10.0   10.0  11.0    10.33   30.17
r4   1.0    1.0   8.5   9.5    4.0    2.33  0.0     1.0     0.83     ;

Table phi(n,n,degree) quadratic fit of the pressure loss (m) div (flow (m^3.h^-1) for each pipe)
                2               1
     s.j1       0.00005425      0.00038190
     j1.j2      0.00027996      0.00149576
     j1.r1      0.00089535      0.00340436
     j1.r4      0.00044768      0.00170218
     j2.r2      0.00223839      0.00851091
     j2.r3      0.00134303      0.00510655;

VARIABLES O                 objectif
          Qk(c,d,t)         debit pompe
          Ql(n,np,t)        debit tuyau entre noeuds n et np
          H(n,t)            hauteur colonne d eau chaque noeud ( = pression(n) .(2*g)^-1 )
          Z(r,t)            pression dans le reservoir
          ;

POSITIVE VARIABLES Qk, H, Ql,Z;
BINARY VARIABLES X(c,d,t);
EQUATIONS OBJ, C4(t,r) , C5(t,c,d), C6(t,n), C7(t), C8(t,n,np) ,C9(c,d,t), C10(r,t);

* conservation debit au niveau reservoir
C4(t,r).. Z(r,t)*surface(r) =E= Z(r,t-1)*surface(r)+SUM(l(n,r),Ql(l,t))-demand(r,t) + vinit(r)$tfirst(t);
* debit pompe limite par un debit max quand elle fonctionne
C5(t,k(c,d)).. Qk(k,t) =L= X(k,t)*Qmax(c);
* Conservation debit a chaque noeud interieur
C6(t,j).. SUM(l(np,j) , Ql(l,t)) =E=  SUM(l(j,np) ,Ql(l,t));
* Conservation du debit au noeud apres les pompes
C7(t)..  SUM(k,Qk(k,t)) =E=  SUM(l('s',np) ,Ql(l,t)) ;
*pertes de charges le long du reseau  CAS GENERAL
C8(t,l(n,np)).. height(n)-height(np) + H(n,t)-H(np,t) =E= phi(l,'1')*Ql(l,t ) + phi(l,'2')*Ql(l,t)*Ql(l,t) ;
* pertes : cas particulier pompes : entre S et j1
C9(k,t).. (H('s',t)-psi('small','0'))*X(k,t) =E= psi('small','2')*Qk(k,t)*Qk(k,t) ;
C10(r,t).. Z(r,t) =L= H(r,t);
* min les couts lie au tarif d electricite pour approvisionner les pompes
OBJ.. O =E= SUM(t,tariff(t)*SUM(k(c,d),gamma(c,'0')*X(k,t)+gamma(c,'1')*Qk(k,t)));

*Bounds
Z.lo(r,t) = vmin(r)/ surface(r);
Z.up(r,t) = vmax(r)/ surface(r);

model M / all /;
solve M using MINLP min O;

display X.l;
display Qk.l;
display Ql.l;
display H.l;     