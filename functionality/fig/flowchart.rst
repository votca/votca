[node distance=.5cm, start chain=going below,]

(init\_global) Global initialization; let 1=(init\_global.north),
2=(init\_global.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode] Initialize global variables (paths to scripts, executables
and user-defined scripts) ;

(init\_step) Iteration initialization; let 1=(init\_step.north),
2=(init\_step.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode] Convert target distribution functions into internal
format, prepare input files, copy data of the previous step;

(init\_sampling) Prepare sampling; let 1=(init\_sampling.north),
2=(init\_sampling.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode] Prepare input files for the external sampling program;

(sampling) Sampling; let 1=(sampling.north), 2=(sampling.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] Canonical ensemble
sampling with molecular dynamics, stochastic dynamics or Monte Carlo
techniques;

(pot\_update) Calculate updates; let 1=(pot\_update.north),
2=(pot\_update.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode, text width=6cm] Analysis of the run. Evaluation of
distribution functions, potential updates :math:`\Delta U^{(n)}` ;

(post\_update) Postprocessing of updates; let 1=(post\_update.north),
2=(post\_update.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode, text width=6cm] Smoothing, extrapolation of potential
updates. Ad-hoc pressure correction. ;

(pot) Update potentials; let 1=(pot.north), 2=(pot.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode, text width=6cm]
:math:`U^{(n+1)} = U^{(n)} + \Delta U^{(n)}` ;

(post\_pot) Postprocessing of potentials; let 1=(post\_pot.north),
2=(post\_pot.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode, text width=6cm] Smoothing, extrapolation of potentials
:math:`U^{(n+1)}` ;

= [diamond, draw = black, very thick, text width=8em, on chain, text
badly centered, node distance=4cm, inner sep=0pt]

(continue) Continue?; let 1=(continue.north), 2=(continue.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode, text width=6cm]
Evaluation of the convergence criterion either for
:math:`\Delta U^{(n)}` or distribution functions. Check the number of
iterations. ;

(end) Finish;

(positive) ;

(continue.west) \|- (positive.east); (positive.north) \|-
(init\_step.west);

(a) at ([xshift=12pt, yshift=-4pt] positive.east) yes; (b) at
([xshift=-4pt] continue.south) no;
