#########################################################
#  Claudia D'Ambrosio, DEIS, University of Bologna
#  15 September 2009
#  nonlinearkp.mod
#########################################################


param pi := 3.14159265358979;

param N; # := 20; # NUMBER OF OBJECTS

set VARS ordered := {1..N};

#option randseed 412836;
#option randseed 412838;
#option randseed 0;

param Umax default 100;
param U {j in VARS}; # default 100;
param a {j in VARS}; # := Uniform(0.1,0.2);
param b {j in VARS}; # := Uniform(0,U);
param c {j in VARS}; # := Uniform(0,U);
param d {j in VARS}; # := Uniform(-U,0);
param weight{VARS};

param C; #default Uniform(1*U,N*U)/1;


var x {j in VARS} >= 0, <= U[j] default Uniform(0,U[j]);
#var x {j in VARS} >= 0, <= U default Uniform(0,U);
var y {VARS};

suffix priority IN, integer, >=0, <=9999;

minimize Total_Profit:
   -sum {j in VARS} y[j];

subject to C_def {j in VARS}:
    y[j] <= c[j]/(1+b[j]*exp(-a[j]*(x[j]+d[j])));

#subject to C_def {j in VARS}:
#   y[j] <= 15*(1/2*sin(2*pi*((x[j]-10)/80)) - cos(pi*((x[j]-10)/80))) + 15*1.3;

subject to KP_constraint:
  sum{j in VARS} weight[j]*x[j] <= C;


