set xlabel "Energy[eV]"
set ylabel "a.u."
plot "vacuum_spec.dat" u 1:2 w impulses lw 5 t "vacuum","static_spec.dat" u 1:2 w impulses lw 5 t "static env."
pause -1
