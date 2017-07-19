%chk=system.chk
%mem=10Gb
%nprocshared=4
#p pop=CHELPG ub3lyp/6-311++G(d,2p) nosymm test polar

TITLE mol 

-1 2
C 30.300016    1.272744    5.521270
H 31.241554    1.478545    6.034844
H 30.490547    0.652427    4.644090
H 29.845289    2.215638    5.209870
H 29.622593    0.750647    6.199925

