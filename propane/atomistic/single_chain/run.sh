#! /bin/bash

grompp -v

mdrun -append -cpi state.cpt &> log_mdrun 

