param a:=10;
set J := {0 to a};
set K := {0 to 1};

var x[J] >=0 <=1;


maximize value: sum <j> in J : floor(random(1,8))*x[j];
subto cons1: forall <l> in K: sum <j> in J : floor(random(0,5))*x[j]<=3.5;

sos s0: type1: x[0]+x[1];
sos s1: type1: x[1]+x[2];
sos s2: type1: x[2]+x[3];
sos s3: type1: x[3]+x[4];
sos s4: type1: x[4]+x[5];
sos s5: type1: x[5]+x[6];
sos s6: type1: x[5]+x[7];
sos s7: type1: x[4]+x[8];
sos s8: type1: x[3]+x[9];
sos s9: type1: x[9]+x[10];