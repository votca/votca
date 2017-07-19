%chk=system.chk
%mem=10Gb
%nprocshared=4
#p pop=CHELPG ub3lyp/6-311++G(d,2p) nosymm test

TITLE mol 

0 1
C 30.300397    1.273753    5.521741
H 31.337631    1.607024    5.782754
H 30.604042    0.445771    4.830929
H 29.841197    2.295753    5.507933
H 29.416733    0.747698    5.966642

