set multiplot layout 2,2

set logscale
set xtics (1,5,10,50,100,500)
set title "convergence"
plot [1:500][0.01:1] "CG-CG.aconv" w l t ""

unset logscale
set xtics autofreq
set title "rdf"
plot [0:0.9][:] "CG-CG.dist.tgt" w l t "target", "CG-CG.dist.new" w l lc 3 t "current"

set title "update"
plot [:][:] "CG-CG.dpot.new" u 1:(4.184*$2) w l t ""

set title "potential"
plot [0:0.9][-2.5:3] "../pot.imc" w l t "imc", "../pot.ibm" w l lc 4 t "ibm" ,"CG-CG.pot.cur" u 1:(4.184*$2) w l lc 3 t "current"

unset multiplot

