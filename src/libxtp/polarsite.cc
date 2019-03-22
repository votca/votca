/*
 *           Copyright 2009-2018 The VOTCA Development Team
 *                      (http://www.votca.org)
 *
 *     Licensed under the Apache License,Version 2.0 (the "License")
 *
 *You may not use this file except in compliance with the License.
 *You may obtain a copy of the License at
 *
 *             http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing,software
 *distributed under the License is distributed on an "AS IS" BASIS,
 *WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,either express or implied.
 *See the License for the specific language governing permissions and
 *limitations under the License.
 *
 */


#include <votca/xtp/polarsite.h>
#include <boost/format.hpp>
#include <fstream>
#include <string>
#include <votca/tools/constants.h>


using namespace std;

namespace votca {
  namespace xtp {

   PolarSite::PolarSite(int id, std::string element, Eigen::Vector3d pos):StaticSite(id,element,pos){
                tools::Elements e;
                double default_pol=std::pow(tools::conv::ang2bohr,3);
                try{
                    default_pol=e.getPolarizability(element)*std::pow(tools::conv::nm2bohr,3);
                }catch(const std::invalid_argument& ){
                   ;
                }
                setPolarisation(default_pol*Eigen::Matrix3d::Identity());
            };


    void PolarSite::Induce(double wSOR) {
      // SUCCESSIVE OVERRELAXATION
      _inducedDipole_old = _inducedDipole; // Remember all previous moments
      _inducedDipole = (1 - wSOR) * _inducedDipole_old - wSOR * _Ps * (_localpermanetField + _localinducedField);
      return;
    }


    std::string PolarSite::writePolarisation()const{
        double conv_pol=std::pow(tools::conv::bohr2ang,3);
        Eigen::MatrixX3d pol=_Ps*conv_pol;
        return (boost::format("     P %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f %6$+1.7f\n")
              % pol(0, 0) % pol(1, 0) % pol(2, 0) % pol(1, 1) % pol(1, 2) % pol(2, 2)).str();
    }


  void PolarSite::SetupCptTable(CptTable& table) const {
      // table.addCol(_Ps[0], "psX", HOFFSET(data, psX));
      // table.addCol(_Ps[1], "psY", HOFFSET(data, psY));
      // table.addCol(_Ps[2], "psZ", HOFFSET(data, psZ));

      table.addCol(_localinducedField[0], "localInducedFieldX", HOFFSET(data, fieldX));
      table.addCol(_localinducedField[1], "localInducedFieldY", HOFFSET(data, fieldY));
      table.addCol(_localinducedField[2], "localInducedFieldZ", HOFFSET(data, fieldZ));

      table.addCol(_inducedDipole[0], "inducedDipoleX", HOFFSET(data, dipoleX));
      table.addCol(_inducedDipole[1], "inducedDipoleY", HOFFSET(data, dipoleY));
      table.addCol(_inducedDipole[2], "inducedDipoleZ", HOFFSET(data, dipoleZ));

      table.addCol(_inducedDipole_old[0], "inducedDipoleXOld", HOFFSET(data, dipoleXOld));
      table.addCol(_inducedDipole_old[1], "inducedDipoleYOld", HOFFSET(data, dipoleYOld));
      table.addCol(_inducedDipole_old[2], "inducedDipoleZOld", HOFFSET(data, dipoleZOld));

      table.addCol(_eigendamp, "eigendamp", HOFFSET(data, eigendamp));
      table.addCol(PhiU, "PhiU", HOFFSET(data, phiU));

  }

  void PolarSite::WriteToCpt(CptTable& table, const std::size_t& idx) const{
      data d;

      d.fieldX     = _localinducedField[0];
      d.fieldY     = _localinducedField[1];
      d.fieldZ     = _localinducedField[2];

      d.dipoleX    = _inducedDipole[0];
      d.dipoleY    = _inducedDipole[1];
      d.dipoleZ    = _inducedDipole[2];

      d.dipoleXOld = _inducedDipole_old[0];
      d.dipoleYOld = _inducedDipole_old[1];
      d.dipoleZOld = _inducedDipole_old[2];

      d.eigendamp = _eigendamp;
      d.phiU = PhiU;

      table.writeToRow(&d, idx);
  }

    void PolarSite::WriteToCpt(const CheckpointWriter& w)const{
        StaticSite::WriteToCpt(w);
        w(_localinducedField,"localinducedField");
        w(_inducedDipole,"inducedDipole");
        w(_inducedDipole_old,"inducedDipole_old");
        w(_eigendamp,"eigendamp");
        w(PhiU,"PhiU");
     }

    void PolarSite::ReadFromCpt(const CheckpointReader& r){
        StaticSite::ReadFromCpt(r);
        r(_localinducedField,"localinducedField");
        r(_inducedDipole,"inducedDipole");
        r(_inducedDipole_old,"inducedDipole_old");
        r(_eigendamp,"eigendamp");
        r(PhiU,"PhiU");
    }

  void PolarSite::ReadFromCpt(CptTable& table, const std::size_t& idx){

      data d;
      table.readFromRow(&d, idx);

      _localinducedField[0] = d.fieldX;
      _localinducedField[1] = d.fieldY;
      _localinducedField[2] = d.fieldZ;

      _inducedDipole[0]     = d.dipoleX;
      _inducedDipole[1]     = d.dipoleY;
      _inducedDipole[2]     = d.dipoleZ;

      _inducedDipole_old[0] = d.dipoleXOld;
      _inducedDipole_old[1] = d.dipoleYOld;
      _inducedDipole_old[2] = d.dipoleZOld;

      _eigendamp            = d.eigendamp;
      PhiU                  = d.phiU;
  }

}
}
