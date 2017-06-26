#!/bin/bash

#calculate bond and angle distributions with csg_boltzmann

cat boltzmann_cmds | csg_boltzmann --top topol.tpr --trj traj.trr --cg hexane.xml

