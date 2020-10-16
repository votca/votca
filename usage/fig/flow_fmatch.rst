[node distance=.5cm, start chain=going below,]

(topol) Reference simulation; let 1=(topol.north), 2=(topol.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] Include forces in
trajectory;

(map) Define mapping scheme; let 1=(map.north), 2=(map.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_dump to list
atoms;

(verify) Verify mapping scheme; let 1=(verify.north), 2=(verify.south)
in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_map to map
Visualize reference + mapped in e.g. VMD;

(input) Setup force-matching options; let 1=(input.north),
2=(input.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode]
Provide correct intervals for distributions (e.g. by csg\_boltzmann,
csg\_stat);

(fmatch) Run force-matching; let 1=(fmatch.north), 2=(fmatch.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_fmatch;

(integrate) Integrate forces to get potential; let 1=(integrate.north),
2=(integrate.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`)
node[tubnode] csg\_call table integrate;
