set I	:= {0 to 8};
param c[I] := <0> 2, <1> 8, <2>4, <3>7, <4>1, <5>9, <6>3,<7>5,<8>7;
var x[I] >=0 <=1;
maximize value: (sum <i> in I : c[i]*x[i]);
subto cons1: 3*x[0]+0*x[1]+9*x[2]+1*x[3]+0*x[4]+3*x[5]+7*x[6]+4*x[7]+6*x[8]<=15;
sos s1: type1: x[0]+x[1];
sos s2: type1: x[0]+x[2];
sos s3: type1: x[1]+x[3];
sos s4: type1: x[1]+x[4];
sos s5: type1: x[3]+x[7];
sos s6: type1: x[3]+x[8];
sos s7: type1: x[2]+x[5];
sos s8: type1: x[2]+x[6];


# /worklocal/scip/scipoptsuite-3.0.1/scip-3.0.1/bin/scip
# read /amd/home/fischer/Desktop/TEst/example.zpl
# optimize
# display solution

# bin/scipmip /amd/home/fischer/Desktop/TEst/example.zpl

