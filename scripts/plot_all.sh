#!/bin/bash

gplot --gnuplot "gnuplot -geometry 400x300+400+0" \
  --term x11 --set logscale \
  --set "xtics (1,5,10,50)" \
  --range [1:55][0.01:1] --title "convergence" \
  convergence.dat w l 

gplot --gnuplot "gnuplot -geometry 400x300+0+0" \
        --term x11 --title "rdf" --range [0:1] \
        -- '"rdf_CG_CG_aim.xvg" w l t "target", "CG-CG.dist.new" w l lc 3 t "current"'

gplot --gnuplot "gnuplot -geometry 400x300+400+350" \
         --term x11 --title "update" CG-CG.dpot.new w l

gplot --gnuplot "gnuplot -geometry 400x300+0+350" \
        --term x11 --range [0:1][-2.5:3] --title "potential" \
       -- '"pot.imc" w l t "imc", "pot.ibm" w l lc 4 t "ibm" ,"CG-CG.pot.new" w l lc 3 t "current"'

