set xlabel "Site energy[eV]"
set ylabel "a.u."

stats "../eanalyze.sitehist_h.out" u 1:2 nooutput

hole_h_mean =STATS_mean_x

stats "../eanalyze.sitehist_e.out" u 1:2 nooutput
hole_e_mean = STATS_mean_x

plot "../eanalyze.sitehist_h.out" u ($1-hole_h_mean):2 w lp t "holes","../eanalyze.sitehist_e.out" u ($1-hole_e_mean):2 w lp t "electrons"
pause -1
