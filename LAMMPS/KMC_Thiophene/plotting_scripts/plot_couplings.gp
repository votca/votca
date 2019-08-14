set xlabel "electronic coupling log[eV**2]"
set ylabel "a.u."


plot "../ianalyze.ihist_h.out" u 1:2 w lp t "holes","../ianalyze.ihist_e.out" u 1:2 w lp t "electrons"
pause -1
