[node distance=.5cm, start chain=going below,]

(ref) Generate target distributions; let 1=(ref.north), 2=(ref.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] Either from
atomistic simulation or experiment;

(topol) Generate coarse-grained topology; let 1=(topol.north),
2=(topol.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode]
either by hand pr csg\_gmxtopol
Cenerate all files to run simulation except for missing potentials;

(input) Generate options file; let 1=(input.north), 2=(input.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] Specify all
interactions that should be iteratively refined;

(ibi) Start iterations; let 1=(ibi.north), 2=(ibi.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_inverse
:math:`<`\ options.xml\ :math:`>`;

(check) Check output; let 1=(check.north), 2=(check.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] Monitor first
couple of iterations.
Many parameters can be tuned on the fly;
