/* 
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#include <string>
#include <boost/filesystem/convenience.hpp> 
#include "dlpolytrajectorywriter.h"

namespace votca { namespace csg {

void DLPOLYTrajectoryWriter::Open(string file, bool bAppend)
{
    if (bAppend)
        throw std::runtime_error("Append to dlpoly trajectory not implemented");

    boost::filesystem::path filepath(file.c_str());
    string filename;
    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        filename="HISTORY_CG";
      } else {
	filename=filepath.parent_path().string() + "/HISTORY_CG";
      }
    } else {
      filename=file;
    }
    _fl.open(filename.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: "+ filename);
}

void DLPOLYTrajectoryWriter::Close()
{
    _fl.close();
}

void DLPOLYTrajectoryWriter::Write(Topology *conf)
{
    static int step = 1;
    int trjkey=0;
    if (conf->HasVel()) trjkey=1; 
    if (conf->HasVel()&& conf->HasForce()) trjkey=2;

    if (step==1) {
      _fl << "From VOTCA with love" << "endl";
      _fl << trjkey << " 2 " << conf->BeadCount() << endl;
    }
    _fl << "timestep 0.0 " << conf->BeadCount() << " " << trjkey << " 0.0 " << conf->getTime();
    matrix m=conf->getBox();
    for (int i=0;i<3;i++) 
      _fl << m[i][1] << " " << m[i][1] << " " << m[i][2] << endl;
    for (int i=0;i<conf->BeadCount(); i++){
      Bead *bead=conf->getBead(i);
      _fl << bead->getName() << " " << i << " " << bead->getM() << " " << bead->getQ() << " 0.0" << endl;
      _fl << bead->getPos().getX() << " " << bead->getPos().getY() << " " << bead->getPos().getZ() << endl;
      if (trjkey>0)
        _fl << bead->getVel().getX() << " " << bead->getVel().getY() << " " << bead->getVel().getZ() << endl;
      if (trjkey>1)
        _fl << bead->getF().getX() << " " << bead->getF().getY() << " " << bead->getF().getZ() << endl;
    }
    step++;
}

}}
