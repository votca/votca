set xlabel "energy [eV]"
set ylabel "a.u."


plot "spectrum.dat" u 1:2 w l t "gaussian","" u 1:4 w l t "lorentzian"
pause -1
