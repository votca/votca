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
        throw std::runtime_error("Error: appending to dlpoly files not implemented");

    boost::filesystem::path filepath(file.c_str());
    string out_name="HISTORY_CGV";

    if ( boost::filesystem::extension(filepath).size() == 0 ) {

      throw std::ios_base::failure("Error on opening transformed dlpoly file '" + file + "' - extension is expected, .dlph or .dlpc");

    } else if( boost::filesystem::extension(filepath)==".dlpc" ) {

      _isConfig=true;
      out_name="CONFIG_CGV";

    } else if( boost::filesystem::extension(filepath)==".dlph" ) {

      _isConfig=false;

    } else {
      throw std::ios_base::failure("Error on opening transformed dlpoly file '" + file + "' - wrong extension, use .dlph or .dlpc");      
    }

    if ( boost::filesystem::basename(filepath).size() == 0 ) {
      if (filepath.parent_path().string().size() == 0) {
        _fname=out_name;
      } else {
        _fname=filepath.parent_path().string() + "/" + out_name;
      }
    } else {
      _fname=file;
    }

    _fl.open(_fname.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on opening transformed dlpoly file '"+ _fname + "'");
}

void DLPOLYTrajectoryWriter::Close()
{
  _fl.close();
}

void DLPOLYTrajectoryWriter::Write(Topology *conf)
{
    static int     step = 1;
    static double dstep = 0.0;
    int    mavecs = 0;
    int    mpbct  = 0;
    double energy = 0.0;

    if( conf->HasForce() && conf->HasVel() ) { mavecs=2; }
    else if               ( conf->HasVel() ) { mavecs=1; }

    if (conf->getBoxType()==BoundaryCondition::typeOrthorhombic) mpbct=2;
    if (conf->getBoxType()==BoundaryCondition::typeTriclinic)    mpbct=3;

    if( _isConfig ) {

      _fl << "From VOTCA with love" << endl;
      _fl << setw(10) << mavecs << setw(10) << mpbct << setw(10) << conf->BeadCount() << setw(20) << energy << endl;
      matrix m=conf->getBox();
      for (int i=0;i<3;i++) 
	_fl << fixed << setprecision(10) << setw(20) << m[i][0] << setw(20) << m[i][1] << setw(20) << m[i][2] << endl;

    } else {

      if (step==1) {
	_fl << "From VOTCA with love" << endl;
	_fl << setw(10) << mavecs << setw(10) << mpbct << setw(10) << conf->BeadCount() << endl;
	dstep = conf->getTime()/(double)(step);
      }

      _fl << "timestep" << setprecision(9) << setw(10) << step << setw(10) << conf->BeadCount() << setw(10) << mavecs << setw(10) << mpbct;
      _fl << setprecision(9) << setw(12) << dstep << setw(12) <<conf->getTime() << endl;

      matrix m=conf->getBox();
      for (int i=0;i<3;i++) 
	_fl << setprecision(12) << setw(20) << m[i][0] << setw(20) << m[i][1] << setw(20) << m[i][2] << endl;

    }

    for (int i=0;i<conf->BeadCount(); i++){
      Bead *bead=conf->getBead(i);

      // AB: DL_POLY needs bead TYPE, not name!

      if( _isConfig) {

	_fl << setw(8) << left << bead->getType()->getName() << right << setw(10) << i+1 << endl;

      } else {

	_fl << setw(8) << left << bead->getType()->getName() << right << setw(10) << i+1;
	_fl << setprecision(6) << setw(12) << bead->getM() << setw(12) << bead->getQ() << setw(12) << "   0.0" << endl;

      }

      // alternative with atom NAME & fixed floating point format (in case the need arises)
      //_fl << setw(8) << left << bead->getName() << right << setw(10) << i+1;
      //_fl << fixed << setprecision(6) << setw(12) << bead->getM() << setw(12) << bead->getQ() << "   0.0" << endl;

      //nm -> Angs
      _fl << resetiosflags(std::ios::fixed) << setprecision(12) << setw(20) << bead->getPos().getX()*10.0;
      _fl << setw(20) << bead->getPos().getY()*10.0 << setw(20) << bead->getPos().getZ()*10.0 << endl;

      if (mavecs>0) {
	if (!bead->HasVel())
	  throw std::ios_base::failure("Error: dlpoly frame is supposed to contain velocities, but bead does not have v-data");

        //nm -> Angs
        _fl << setprecision(12) << setw(20) << bead->getVel().getX()*10.0 << setw(20);
        _fl << bead->getVel().getY()*10.0 << setw(20) << bead->getVel().getZ()*10.0 << endl;

	if (mavecs>1) {
	  if (!bead->HasF())
	    throw std::ios_base::failure("Error: dlpoly frame is supposed to contain forces, but bead does not have f-data");

          //nm -> Angs
	  _fl << setprecision(12) << setw(20) << bead->getF().getX()*10.0 << setw(20);
	  _fl << bead->getF().getY()*10.0 << setw(20) << bead->getF().getZ()*10.0 << endl;
	}
      }
    }
    step++;
}

}}

