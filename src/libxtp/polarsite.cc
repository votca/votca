/* 
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
/// For earlier commit history see ctp commit 77795ea591b29e664153f9404c8655ba28dc14e9

#include <votca/xtp/polarsite.h>
#include <boost/math/special_functions/round.hpp>
#include <boost/format.hpp>
#include <fstream>
#include <string>


using namespace std;

namespace votca { namespace xtp {


 

void PolarSite::Rotate(const tools::matrix &rot, const tools::vec &refPos) {

    tools::vec dir = _pos - refPos;
    dir = rot * dir;
    _pos = refPos + dir;
       
    tools::matrix R_T = tools::matrix(rot.getRow(0),rot.getRow(1),rot.getRow(2));
        
        // Transform polarizability tensor into global frame
        tools::matrix P_Global = R * _Ps * R_T;
        _Ps = P_Global;
        
        // Any multipoles for this charge state available?
        if (_multipoles.size() < 1) { continue; }

        // Transform dipole moment into global frame
        if (_multipoles.size() > 1) {

            double Qz = _Qs[state+1][1];
            double Qx = _Qs[state+1][2];
            double Qy = _Qs[state+1][3];

            vec d = vec(Qx, Qy, Qz);
            d = R * d;

            _Qs[state+1][1] = d.getZ();
            _Qs[state+1][2] = d.getX();
            _Qs[state+1][3] = d.getY();
        }

        // Transform quadrupole moment into global frame
        if (_multipoles.size() > 4) {

            double Qzz =      _Qs[state+1][4];
            double Qxx = -0.5*_Qs[state+1][4] + 0.5*sqrt(3)*_Qs[state+1][7];
            double Qyy = -0.5*_Qs[state+1][4] - 0.5*sqrt(3)*_Qs[state+1][7];

            double Qxy =  0.5*sqrt(3)*_Qs[state+1][8];
            double Qxz =  0.5*sqrt(3)*_Qs[state+1][5];
            double Qyz =  0.5*sqrt(3)*_Qs[state+1][6];

            matrix Q = matrix(vec(Qxx,Qxy,Qxz),
                              vec(Qxy,Qyy,Qyz),
                              vec(Qxz,Qyz,Qzz));

            matrix Q_Global  = R * Q * R_T;
            
      
            _Qs[state+1][4] =               Q_Global.get(2,2);  // Q20
            _Qs[state+1][5] = 2 / sqrt(3) * Q_Global.get(0,2);  // Q21c
            _Qs[state+1][6] = 2 / sqrt(3) * Q_Global.get(1,2);  // Q21s
            _Qs[state+1][7] = 1 / sqrt(3) *(Q_Global.get(0,0)
                                          - Q_Global.get(1,1)); // Q22c
            _Qs[state+1][8] = 2 / sqrt(3) * Q_Global.get(0,1);  // Q22s
        }

    }

return;
}



void PolarSite::Translate(const tools::vec &shift) {
    _pos += shift;
    return;
}

void PolarSite::Induce(double wSOR) {
    // SUCCESSIVE OVERRELAXATION
    _inducedDipole_old=_inducedDipole;// Remember all previous moments
    
    _inducedDipole=(1-wSOR)*_inducedDipole_old-_Ps*(_localpermanetField+_localinducedField);
    return;  
}

void PolarSite::WriteMpsLine(std::ostream &out, string unit = "angstrom") {
    
    // Set conversion factor for higher-rank moments (e*nm**k to e*a0**k)
    double conv_dpl = 1./0.0529189379;
    double conv_qdr = conv_dpl*conv_dpl;
    // Set conversion factor for polarizabilities (nm**3 to A**3)
    double conv_pol = 1000;    
    // Set conversion factor for positions (nm to ??)
    double conv_pos = 1.;
    if (unit == "angstrom") {
        conv_pos = 10.;
    }
    else if (unit == "nanometer") {
        conv_pos = 1.;
    }
    else assert(false); // Units error
    
    out << (boost::format(" %1$2s %2$+1.7f %3$+1.7f %4$+1.7f Rank %5$d\n") 
            % _name % (_pos.getX()*conv_pos)
            % (_pos.getY()*conv_pos) % (_pos.getZ()*conv_pos)
            % _rank);
    // Charged
    out << (boost::format("    %1$+1.7f\n") % Q00);
    if (_rank > 0) {
        // Dipole z x y
        out << (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f\n") 
            % (_multipoles(3)*conv_dpl) % (_multipoles(1)*conv_dpl) % (_multipoles(2)*conv_dpl));
        if (_rank > 1) {
            // Quadrupole 20 21c 21s 22c 22s
            out << (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f\n") 
                % (_multipoles(4)*conv_qdr) % (_multipoles(5)*conv_qdr) % (_multipoles(6)*conv_qdr) 
                % (_multipoles(7)*conv_qdr) % (_multipoles(8)*conv_qdr));
        }
    }
    // Polarizability
    out << (boost::format("     P %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f %6$+1.7f \n") 
        % (Pxx*conv_pol) % (Pxy*conv_pol) % (Pxz*conv_pol) 
        % (Pyy*conv_pol) % (Pyz*conv_pol) % (Pzz*conv_pol));
    
}




}}
