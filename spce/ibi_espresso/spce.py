#!/usr/bin/env python3

# DON'T RUN THIS MANUALLY
# RUN "csg_inverse --options settings.xml" instead

import time
import numpy as np

import espressomd
from espressomd.io.writer import h5md

def readgrofile(filename):
    atompos=[]
    box=[]
    atomnumber=0
    with open(filename,"r") as f:
        comment=f.readline()
        atomnumber=int(f.readline())
        for i in range(atomnumber):
            line=f.readline()
            atompos.append(line.split()[3:])
        box=(f.readline().split())
    return atomnumber,np.array(box,dtype=float),np.array(atompos,dtype=float)
    

def write_data(file_name, time, energy, n_part):
    with open(file_name, 'a') as f:
        e = energy['total'] / n_part
        e_pot = (energy['total'] - energy['kinetic']) / n_part
        e_kin = e - e_pot
        print(" {} {} {} {}".format(time, e, e_pot, e_kin), file=f)


def calc_temperature(system):
    E = system.analysis.energy()['total']
    n_part = len(system.part)

    return 2. / 3. * E / n_part

# check for necessary feature
required_features = ["MASS", "TABULATED"]
espressomd.assert_features(required_features)

print("Program Information:\n{}\n".format(espressomd.features()))

# set system properties:
skin = 0.5
time_step = 0.002
friction = 5.0
int_steps = 900
eq_steps = 100
steps_per_int = 100

atomnumber,box,atompos=readgrofile('spce.gro')
nm2Angstroem=10
box_length = box*nm2Angstroem
positions = atompos*nm2Angstroem
mass=18
masses=mass*np.ones(atomnumber)

# Setup Espresso with particle from spce.gro
system = espressomd.System(box_l=box_length, time_step=time_step)
system.cell_system.skin = skin
system.set_random_state_PRNG()
system.thermostat.set_langevin(kT=1., gamma=1.)
system.part.add(pos=positions, mass=masses)

print("number of particles: ", len(system.part))
print("box size: ", box_length)
print("initial temperature: ", calc_temperature(system))

with open('CG_CG.tab', 'r') as f:
    data = np.loadtxt(f)
    r = data[:, 0]
    f = data[:, 1]
    p = data[:, 2]

    min_r = np.min(r)
    max_r = np.max(r)

    system.non_bonded_inter[0, 0].tabulated.set_params(min=min_r, max=max_r,
                                                       energy=p,
                                                       force=f)

# Warmup loop
for i in range(eq_steps):
    print("warmup step", i)
    system.integrator.run(steps_per_int)

print("Running at temperature T={:.2f}".format(calc_temperature(system)))

h5_file = h5md.H5md(filename="traj.h5", write_pos=True, write_vel=True,
                    write_force=True, write_species=False, write_mass=True, write_ordered=False)

starttime = time.time()

# Integration loop
for i in range(int_steps):
    system.integrator.run(steps_per_int)

    E = system.analysis.energy()
    e_pot = (E['total'] - E['kinetic']) / len(system.part)
    print("time: {:.3f} potential energy: {:.2f}".format(i * time_step, e_pot))

    h5_file.write()
    write_data("energy.dat", i * time_step,
               system.analysis.energy(), len(system.part))

endtime = time.time()
h5_file.close()

print("time per MD step: ", 1000. * (endtime - starttime) / int_steps, "ms")
print("final temperature: ", calc_temperature(system))
print("Finished.")
