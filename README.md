# BD
solve covering problem by BD in GAMS
********************************************************
************** BD  *************************************
********************************************************

Sets
I /i1*i20/
J /j1*j100/
alias (i,ip)
;

*********** Data
Parameters
B
cap(i)
c(i)
dem(j)
w(j)
dis(i,j)
r
;

B        = 1000;
cap(i)   = uniform(100,200);
c(i)     = uniform(100,200);
dem(j)   = uniform(20,50);
w(j)     = uniform(0,1);
dis(i,j) = uniform(5,50);
r        = 15;

Parameter
dis_check
;

dis_check(i,j)=0;
dis_check(i,j)$(dis(i,j)/r<1)=1;

Display


B
cap
c
dem
w
dis
r

********************************************************************************
* Step 0: Initialization of BD**************************************************

Set iter /1*20/

Parameters
MaxRE   /1E-3/
LB
UB
;
LB=0.0000001;
UB=inf;

********************************************************************************

********************************************************************************
* Step 1: Subproblem (SP) of BD*************************************************

Free variables
Z_SP    ;

Positive variable
y(i,j) ;

Parameter
x_fix(i) ;

Equations
obj_SP       'obective function of SP'
Cons1_SP      'Constraints of SP'
Cons2_SP
Cons3_SP
Cons4_SP
;

obj_SP..         Z_SP =e= sum({i,j},w(j)*y(i,j));

Cons1_SP       ..  sum({i},c(i)*x_fix(i)) =l= B ;


Cons2_SP(i)    ..  sum({j},y(i,j)) =l= cap(i)*x_fix(i);


Cons3_SP(i,j)  ..  y(i,j) =l= dem(j)*dis_check(i,j);


Cons4_SP(j)    ..  sum({i},y(i,j)) =l= dem(j);
 


Model SP
/
obj_SP
Cons1_SP
Cons2_SP
Cons3_SP
Cons4_SP
/
;

********************************************************************************
********************************************************************************
********************************************************************************
* Step 2: Dual Subproblem (DSP) of BD*******************************************

Free variables
Z_DSP    ;

Positive variables
alfa(i)
beta(j)
gama(i,j) ;


Equations
obj_DSP       'obective function of DSP'
Cons_DSP      'Constraints of DSP'
;

obj_DSP..         Z_DSP =e=  sum({i},cap(i)*x_fix(i)*alfa(i))+sum({j},dem(j)*beta(j))+sum((i,j),dem(j)*gama(i,j));

Cons_DSP(i,j)..   alfa(i) + beta(j) + gama(i,j )=g= w(j);

Model DSP
/
obj_DSP
Cons_DSP
/
;
********************************************************************************
********************************************************************************
********************************************************************************
* Step 3: Homogeneous/Modified Dual Subproblem (MDSP) of BD*********************

Free variables
Z_MDSP    ;

Equations
obj_MDSP       'obective function of MDSP'
Cons_MDSP      'Constraints of MDSP'
Bounded_Cons1
Bounded_Cons2
Bounded_Cons3
;

obj_MDSP..         Z_MDSP =e=  sum({i},cap(i)*x_fix(i)*alfa(i))+sum({j},dem(j)*beta(j))+sum((i,j),dem(j)*gama(i,j));
Cons_MDSP(i,j)..     alfa(i) + beta(j) + gama(i,j) =g= 0 ;
Bounded_Cons1..    sum(i,alfa(i)) =g= card(i);
Bounded_Cons2..    sum(j,beta(j)) =g= card(j);
Bounded_Cons3..    sum((i,j),gama(i,j)) =g= card(i)+card(j);
Model MDSP
/
obj_MDSP
Cons_MDSP
Bounded_Cons1
Bounded_Cons2
Bounded_Cons3
/
;
********************************************************************************
********************************************************************************
********************************************************************************
* Step 4: Relaxed Master Problem (RMP) of BD************************************

Free variables
Z_RMP    ;

Binary variable
x(i)
;

Positive Variable Say;

Sets
OC(iter)
FC(iter)
;

OC(iter)= NO;
FC(iter)= NO;

Parameters
alfa_fix(i,iter)
beta_fix(j,iter)
gama_fix(i,j,iter)
;



Equations
obj_RMP       'obective function of MP'


OptimalityCut(iter)
FeasibilityCut(iter)

;

obj_RMP..                Z_RMP =e= Say;


OptimalityCut(OC)..      Say  =l= sum({i},cap(i)*x(i)*alfa_fix(i,oc))+sum({j},dem(j)*beta_fix(j,oc))+sum((i,j),dem(j)*gama_fix(i,j,oc));

FeasibilityCut(FC) ..   sum({i},cap(i)*x(i)*alfa_fix(i,fc))+sum({j},dem(j)*beta_fix(j,fc))+sum((i,j),dem(j)*gama_fix(i,j,fc)) =g= 0;


Model RMP
/
obj_RMP

OptimalityCut
FeasibilityCut
/
;
********************************************************************************
* Step 5: Main Loop for implentation of BD**************************************


Parameters
Result(iter,*)
Converged /NO/
Iteration
Gap
x_Feasibility
;

x_fix(i)=0;

Options
LP  =  CPLEX
MIP = CPLEX
OPTCR = 0
RESLIM = 300
;

Loop(iter$(NOT(Converged)),
***** Solve DSP or MDSP to find u and update LB
Solve DSP us LP Min Z_DSP;

*Infeasibility of Main OP
Abort$(DSP.ModelStat = 2) "Your OP Model is not fasible"

* Bounded Situtation
if( DSP.ModelStat <> 3,

x_Feasibility = YES;
OC(iter)=YES;
LB= Z_DSP.l;
Result(iter,'LB')=LB;

ELSE
* Unbounded Situtation
x_Feasibility = NO;
FC(iter)=YES;
Solve MDSP us LP Min Z_MDSP;

)
* end of If
;

Result(iter,'Feasible')=x_Feasibility;
alfa_fix(i,iter)=alfa.l(i);
beta_fix(j,iter)=beta.l(j);
gama_fix(i,j,iter)=gama.l(i,j);
************************************************
*************** Solve RMP to  fine new y and y and update UB

Solve RMP us MIP Max Z_RMP ;

Abort$(RMP.ModelStat = 2) "Your OP Model is not fasible" ;

UB=Z_RMP.l;
Result(iter,'UB')=UB;

x_fix(i)=x.l(i);
******************************************


* Stop Criteria

Gap = abs((UB - LB)/ LB );

Result(iter,'Gap')=Gap;

IF(Gap <= MaxRE,
Converged = YES;
)
;

Iteration=ord(iter);
Display
"-----------Iteration-----------"
Iteration
x_Feasibility
LB
UB
Gap
"-----------Varable-------------"
"x"
x_fix
"y"
*Cons_DSP.i,j

)
;
*End of BD Main Loop


Display
Result

