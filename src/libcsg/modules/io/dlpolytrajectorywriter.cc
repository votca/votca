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

#include <iomanip>
#include <string>
#include <boost/filesystem/convenience.hpp> 
#include "dlpolytrajectorywriter.h"

namespace votca { namespace csg {

void DLPOLYTrajectoryWriter::Open(string file, bool bAppend)
{
    if (bAppend)
        throw std::runtime_error("Append to dlpoly HISTORY (trajectory) not implemented");

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
      throw std::ios_base::failure("Error on open dlpoly file: '" + filename + "'");
}

void DLPOLYTrajectoryWriter::Close()
{
  _fl.close();
}

void DLPOLYTrajectoryWriter::Write(Topology *conf)
{
    static int step = 1;
    int mavecs=0;
    int mpbct=0;

    if      (conf->HasForce()) { mavecs=2; }
    else if (conf->HasVel())   { mavecs=1; }

    if (conf->getBoxType()==BoundaryCondition::typeOrthorhombic) mavecs=2;
    if (conf->getBoxType()==BoundaryCondition::typeTriclinic)    mavecs=3;

    if (step==1) {
      _fl << "From VOTCA with love" << endl;
      _fl << setw(10) << mavecs << setw(10) << mpbct << setw(10) << conf->BeadCount() << endl;
    }

    _fl << "timestep" << setprecision(9) << setw(10) << step << setw(10) << conf->BeadCount() << setw(10) << mavecs << setw(10) << mpbct;
    _fl << setprecision(16) << setw(20) << conf->getTime()/(double)(step) << setw(20) <<conf->getTime() << endl;

    matrix m=conf->getBox();
    for (int i=0;i<3;i++) 
      _fl << fixed << setprecision(10) << setw(20) << m[i][0] << setw(20) << m[i][1] << setw(20) << m[i][2] << endl;

    for (int i=0;i<conf->BeadCount(); i++){
      Bead *bead=conf->getBead(i);
      _fl << bead->getName() << setprecision(9) << setw(10) << i << setw(12) << bead->getM() << setw(12) << bead->getQ() << "   0.0" << endl;
      _fl << setprecision(12) << setw(16) << bead->getPos().getX() << setw(16) << bead->getPos().getY() << setw(16) << bead->getPos().getZ() << endl;
      if (mavecs>0)
        _fl << setprecision(12) << setw(16) << bead->getVel().getX() << setw(16) << bead->getVel().getY() << setw(16) << bead->getVel().getZ() << endl;
      if (mavecs>1)
        _fl << setprecision(12) << setw(16) << bead->getF().getX() << setw(16) << bead->getF().getY() << setw(16) << bead->getF().getZ() << endl;
    }
    step++;
}

}}

