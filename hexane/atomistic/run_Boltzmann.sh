#!/bin/bash -e

#calculate bond and angle distributions with csg_boltzmann

csg_boltzmann --top topol.tpr --trj traj.trr --cg hexane.xml < boltzmann_cmds

