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

   PolarSite::PolarSite(int id, std::string element, Eigen::Vector3d pos)
            : _id(id), _element(element),_pos(pos),_rank(0),_multipole(Eigen::VectorXd::Zero(1)),
           _localpermanetField(Eigen::Vector3d::Zero()),
            _localinducedField(Eigen::Vector3d::Zero()),_eigendamp(0.0),
            PhiP(0.0),PhiU(0.0){
               
                tools::Elements e;
                double default_pol=std::pow(tools::conv::ang2bohr,3);
                try{
                    default_pol=e.getPolarizability(element)*std::pow(tools::conv::nm2bohr,3);
                }catch(const std::invalid_argument& ){
                     std::cout << std::endl << "WARNING No default polarizability given for "
                    << "polar site type '" << element << "'. Defaulting to 1 A**3. "
                    << std::flush;
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

  }
}
