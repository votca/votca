/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <stdio.h>
#include <string>
#include "lammpsdumpwriter.h"

namespace votca { namespace csg {

using namespace std;

void LAMMPSDumpWriter::Open(std::string file, bool bAppend)
{
    _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void LAMMPSDumpWriter::Close()
{
    fclose(_out);
}

void LAMMPSDumpWriter::Write(Topology *conf)
{
    Topology *top = conf;
    votca::tools::matrix box = conf->getBox();
    fprintf(_out, "ITEM: TIMESTEP\n%i\n", top->getStep());
    fprintf(_out, "ITEM: NUMBER OF ATOMS\n%i\n", (int)top->Beads().size());
    fprintf(_out, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(_out, "0 %f\n0 %f\n0 %f\n",box[0][0],box[1][1],box[2][2]);

    fprintf(_out, "ITEM: ATOMS id type x y z");
    bool v = top->HasVel();
    if(v) {
      fprintf(_out, " vx vy vz");
    }
    bool f = top->HasForce();
    if(f) {
      fprintf(_out, " fx fy fz");
    }
    fprintf(_out, "\n");

    for(BeadContainer::iterator iter=conf->Beads().begin(); iter!=conf->Beads().end(); ++iter) {
        Bead *bi = *iter;
        fprintf(_out,"%i %i", bi->getId()+1, bi->getType()->getId());
        fprintf(_out," %f %f %f",bi->getPos().getX(), bi->getPos().getY(), bi->getPos().getZ());
        if(v) {
            fprintf(_out, " %f %f %f", bi->getVel().getX(), bi->getVel().getY(), bi->getVel().getZ());
        }
        if(f) {
            fprintf(_out, " %f %f %f", bi->getF().getX(), bi->getF().getY(), bi->getF().getZ());
        }
        fprintf(_out, "\n");
    }
    fflush(_out);
}

}}
