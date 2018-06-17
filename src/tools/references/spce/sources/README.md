Gromacs input data, together with spc216.gro from Gromacs distributions.
```
cp /usr/share/gromacs/top/spc216.gro .
gmx grompp -c spc216.gro
gmx mdrun
csg_map --top topol.tpr --trj traj.trr --first-frame=2 --no-map --out frame.dump --force
```
