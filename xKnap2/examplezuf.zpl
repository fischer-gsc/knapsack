param a:=30;
param d:=5;

set J	:= {0 to a};
set K	:= {1 to a};
set L	:= {0 to d};
var x[J] >=0 <=1;

#maximize value: sum <j> in J : floor(floor(random(0,d))*random(0,a))*x[j];
#subto cons1: sum <j> in J : random(0,a)*x[j]<=random(1*a,10*a);

maximize value: sum <j> in J : floor(random(0,a))*x[j];
subto cons1: forall <l> in L: sum <j> in J : random(0,a)*x[j]<=90.523;
sos sos1_cons1: forall <k> in K: type1: 1*x[floor(random(0,k-1+0.00000000001))] + 1*x[k-1];


# /worklocal/scip/scipoptsuite-3.0.1/scip-3.0.1/bin/scip

# read /amd/home/fischer/Desktop/TEst/example.zpl

# optimize

# display solution

# bin/scipmip /amd/home/fischer/Desktop/TEst/example.zpl

