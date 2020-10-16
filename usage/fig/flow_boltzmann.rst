[node distance=.5cm, start chain=going below,]

(topol) Prepare atomistic topology; let 1=(topol.north), 2=(topol.south)
in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] ;

(map) Define mapping scheme; let 1=(map.north), 2=(map.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_dump to list
atoms;

(verify) Verify mapping scheme; let 1=(verify.north), 2=(verify.south)
in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_map to map
Visualize reference + mapped in e.g. VMD;

(excl) Create exclusion list; let 1=(excl.north), 2=(excl.south) in
(:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] csg\_boltzmann
–excl;

(traj) Generate reference trajectory; let 1=(traj.north), 2=(traj.south)
in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] ;

(pot) csg\_boltzmann to get distributions/potentials; let 1=(pot.north),
2=(pot.south) in (:math:`(2, \y1)`) – (:math:`(2, \y2)`) node[tubnode] ;
