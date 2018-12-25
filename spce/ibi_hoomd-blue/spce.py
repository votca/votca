#!/usr/bin/env python2

# DON'T RUN THIS MANUALLY
# RUN "csg_inverse --options settings.xml" instead

from hoomd_script import *

'''
2180 Particles
Box 4.031^3
langevin thermostat T=2.5 gamma=5
dt=0.002
cutoff=0.9
'''

def create_system_from_gro(gro_file):
    """ Create system using position from a GROMACS gro file. """

    f = open(gro_file)
    f.readline() # skip comment line
    n_particles = int(f.readline())
    
    # store coordinates and velocities
    x, y, z = [], [], []
    for i in range(n_particles):
        s = f.readline()[20:69]
        # coordinates
        x.append(float(s[0:8]))
        y.append(float(s[8:16]))
        z.append(float(s[16:24]))
    
    # store box size
    Lx, Ly, Lz = map(float, f.readline().split()) # read last line, convert to float
    f.close()

    system = init.create_empty(N=n_particles, box=data.boxdim(Lx=Lx,Ly=Ly,Lz=Lz), particle_types=['CG'])
    for i, p in enumerate(system.particles):
        p.position = ( x[i], y[i], z[i] )  
        p.type = 'CG'

    return(system)


def write_gro(gro_file):
    """ Write positions to GROMACS gro file. """
    
    f = open(gro_file,'a')
    st = "gro write by hoomd-blue step=%i\n"%get_step()
    f.write(st)
    f.write(' '+str(len(system.particles))+'\n')
    for i, p in enumerate(system.particles):
        #st = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"%(i,'NON',p.type,p.postion,p.velocity)
        st = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"%(i+1,'SOL',p.type,i+1,
            p.position[0],p.position[1],p.position[2],p.velocity[0],p.velocity[1],p.velocity[2])
        f.write(st)
    st = "%10.5f%10.5f%10.5f\n"%(system.box.Lx,system.box.Ly,system.box.Lz)
    f.write(st)
    f.close()


system = create_system_from_gro('conf.gro')
tab = pair.table(width=351)
tab.set_from_file('CG', 'CG', filename="CG_CG.tab")

all = group.all();
integrate.mode_standard(dt=0.002)
bd = integrate.bdnvt(group=all, T=2.5)
bd.set_gamma('CG', gamma=5.0)

for step in range(100):
  run(100)

for step in range(900):
  run(100)
  write_gro('traj.gro')

print "finished"
