%chk=system.chk
%mem=10Gb
%nprocshared=4
#opt pop=CHELPG ub3lyp/6-311++G(d,2p) nosymm test

TITLE mol 

1 2
C         30.30000        1.27000        5.52000
H         31.24000        1.48000        6.03000
H         30.49000        0.66000        4.65000
H         29.85000        2.21000        5.21000
H         29.62000        0.75000        6.20000

