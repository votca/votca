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
      table.addCol(_id, "index", HOFFSET(data, id));
      table.addCol(_element, "type", HOFFSET(data, element));

      table.addCol(_pos[0], "posX", HOFFSET(data, posX));
      table.addCol(_pos[1], "posY", HOFFSET(data, posY));
      table.addCol(_pos[2], "posZ", HOFFSET(data, posZ));

      table.addCol(_rank, "rank", HOFFSET(data, rank));

      table.addCol(_multipole[0], "multipoleQ00",  HOFFSET(data, multipoleQ00));
      table.addCol(_multipole[1], "multipoleQ11c", HOFFSET(data, multipoleQ11c));
      table.addCol(_multipole[2], "multipoleQ11s", HOFFSET(data, multipoleQ11s));
      table.addCol(_multipole[3], "multipoleQ10",  HOFFSET(data, multipoleQ10));
      table.addCol(_multipole[4], "multipoleQ20",  HOFFSET(data, multipoleQ20));
      table.addCol(_multipole[5], "multipoleQ21c", HOFFSET(data, multipoleQ21c));
      table.addCol(_multipole[6], "multipoleQ21s", HOFFSET(data, multipoleQ21s));
      table.addCol(_multipole[7], "multipoleQ22c", HOFFSET(data, multipoleQ22c));
      table.addCol(_multipole[8], "multipoleQ22s", HOFFSET(data, multipoleQ22s));

      table.addCol(_Ps(0,0), "pxx",  HOFFSET(data, pxx));
      table.addCol(_Ps(0,1), "pxy",  HOFFSET(data, pxy));
      table.addCol(_Ps(0,2), "pxz",  HOFFSET(data, pxz));
      table.addCol(_Ps(1,1), "pyy",  HOFFSET(data, pyy));
      table.addCol(_Ps(1,2), "pyz",  HOFFSET(data, pyz));
      table.addCol(_Ps(2,2), "pzz",  HOFFSET(data, pzz));

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

      d.id = _id;
      d.element = const_cast<char*>(_element.c_str());
      d.posX = _pos[0];
      d.posY = _pos[1];
      d.posZ = _pos[2];

      d.rank = _rank;

      d.multipoleQ00  = _multipole[0];
      d.multipoleQ11c = _multipole[1];
      d.multipoleQ11s = _multipole[2];
      d.multipoleQ10  = _multipole[3];
      d.multipoleQ20  = _multipole[4];
      d.multipoleQ21c = _multipole[5];
      d.multipoleQ21s = _multipole[6];
      d.multipoleQ22c = _multipole[7];
      d.multipoleQ22s = _multipole[8];

    d.pxx=_Ps(0,0);
    d.pxy=_Ps(0,1);
    d.pxz=_Ps(0,2);
    d.pyy=_Ps(1,1);
    d.pyz=_Ps(1,2);
    d.pzz=_Ps(2,2);

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

  void PolarSite::WriteData(data& d) const{
      d.id = _id;
      d.element = const_cast<char*>(_element.c_str());
      d.posX = _pos[0];
      d.posY = _pos[1];
      d.posZ = _pos[2];

      d.rank = _rank;

      d.multipoleQ00  = _multipole[0];
      d.multipoleQ11c = _multipole[1];
      d.multipoleQ11s = _multipole[2];
      d.multipoleQ10  = _multipole[3];
      d.multipoleQ20  = _multipole[4];
      d.multipoleQ21c = _multipole[5];
      d.multipoleQ21s = _multipole[6];
      d.multipoleQ22c = _multipole[7];
      d.multipoleQ22s = _multipole[8];

      d.pxx=_Ps(0,0);
      d.pxy=_Ps(0,1);
      d.pxz=_Ps(0,2);
      d.pyy=_Ps(1,1);
      d.pyz=_Ps(1,2);
      d.pzz=_Ps(2,2);

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
  }

  void PolarSite::ReadFromCpt(CptTable& table, const std::size_t& idx){

      data d;
      table.readFromRow(&d, idx);
      _id           = d.id;
      _element      = std::string(d.element);
      free(d.element);
      _pos[0]       = d.posX;
      _pos[1]       = d.posY;
      _pos[2]       = d.posZ;

      _rank         = d.rank;

      _multipole[0] = d.multipoleQ00;
      _multipole[1] = d.multipoleQ11c;
      _multipole[2] = d.multipoleQ11s;
      _multipole[3] = d.multipoleQ10;
      _multipole[4] = d.multipoleQ20;
      _multipole[5] = d.multipoleQ21c;
      _multipole[6] = d.multipoleQ21s;
      _multipole[7] = d.multipoleQ22c;
      _multipole[8] = d.multipoleQ22s;

      _localinducedField[0] = d.fieldX;
      _localinducedField[1] = d.fieldY;
      _localinducedField[2] = d.fieldZ;

      _inducedDipole[0]     = d.dipoleX;
      _inducedDipole[1]     = d.dipoleY;
      _inducedDipole[2]     = d.dipoleZ;

      _inducedDipole_old[0] = d.dipoleXOld;
      _inducedDipole_old[1] = d.dipoleYOld;
      _inducedDipole_old[2] = d.dipoleZOld;

    _Ps(0,0)=d.pxx;
    _Ps(0,1)=d.pxy;
    _Ps(1,0)=d.pxy;
    _Ps(0,2)=d.pxz;
    _Ps(2,0)=d.pxz;
    _Ps(1,1)=d.pyy;
    _Ps(1,2)=d.pyz;
    _Ps(2,1)=d.pyz;
    _Ps(2,2)=d.pzz;

      _eigendamp            = d.eigendamp;
      PhiU                  = d.phiU;
  }

  void PolarSite::ReadData(data& d){
      _id           = d.id;
      _element      = std::string(d.element);
      free(d.element);
      _pos[0]       = d.posX;
      _pos[1]       = d.posY;
      _pos[2]       = d.posZ;

      _rank         = d.rank;

      _multipole[0] = d.multipoleQ00;
      _multipole[1] = d.multipoleQ11c;
      _multipole[2] = d.multipoleQ11s;
      _multipole[3] = d.multipoleQ10;
      _multipole[4] = d.multipoleQ20;
      _multipole[5] = d.multipoleQ21c;
      _multipole[6] = d.multipoleQ21s;
      _multipole[7] = d.multipoleQ22c;
      _multipole[8] = d.multipoleQ22s;

      _localinducedField[0] = d.fieldX;
      _localinducedField[1] = d.fieldY;
      _localinducedField[2] = d.fieldZ;

      _inducedDipole[0]     = d.dipoleX;
      _inducedDipole[1]     = d.dipoleY;
      _inducedDipole[2]     = d.dipoleZ;

      _inducedDipole_old[0] = d.dipoleXOld;
      _inducedDipole_old[1] = d.dipoleYOld;
      _inducedDipole_old[2] = d.dipoleZOld;

      _Ps(0,0)=d.pxx;
      _Ps(0,1)=d.pxy;
      _Ps(1,0)=d.pxy;
      _Ps(0,2)=d.pxz;
      _Ps(2,0)=d.pxz;
      _Ps(1,1)=d.pyy;
      _Ps(1,2)=d.pyz;
      _Ps(2,1)=d.pyz;
      _Ps(2,2)=d.pzz;

      _eigendamp            = d.eigendamp;
      PhiU                  = d.phiU;
  }

}
}
