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

// Get Rank: We have at the moment just three cases. Point charges, dipoles and Quadrupoles
    void PolarSite::calcRank(){
        int mplen = _multipoles.size();
        if(mplen == 1){_rank = 0;}
        else if(mplen == 4){_rank = 1;}
        else if(mplen == 9){_rank = 2;}
        else{
            throw std::runtime_error("PolarSite. multipoles size is not 1,4 or 9.");
        }
        return;
    }
    
// This function put in the conventional order the dipoles
    // From Spherical (Q_10=mu_z , Q_11c=mu_x, Q_11s=mu_y)-->(Q_11c,Q_11s,Q_10)
Eigen::Vector3d PolarSite::getCartesianDipoles(){
    Eigen::Vector3d cartesiandipole;
    if(_multipoles.size() > 1){
    Eigen::VectorXd  MP = _multipoles;
    cartesiandipole(0)=MP(2);
    cartesiandipole(1)=MP(3);
    cartesiandipole(2)=MP(1);
    return cartesiandipole;
    }
}    
    //Transform Multipoles from Spherical to Cartesian
    // spherical_multipoles Q = ( Q00,Q10,Q11c,Q11s,Q20,Q21c,Q21s,Q22c,Q22s )
    // We are trasforming here just quadrupoles
 Eigen::Matrix3d PolarSite::getCartesianMultipoles() {
     Eigen::VectorXd MP = _multipoles;
     if( _rank > 1 ){
        Eigen::Matrix3d theta;
        theta(0,0) = 0.5 * (-MP(4)+std::sqrt(3.) * MP(7)); // theta_xx
        theta(1,1) = 0.5 * (-MP(4)+std::sqrt(3.) * (-MP(7))); // theta_yy
        theta(2,2) = MP(4); // theta_zz
        theta(0,1) = theta(1,0) = 0.5 * std::sqrt(3.) * MP(6); // theta_xy = theta_yx
        theta(0,2) = theta(2,0) = 0.5 * std::sqrt(3.) * MP(5); // theta_xz = theta_zx
        theta(1,2) = theta(2,1) = 0.5 * std::sqrt(3.) * MP(6); //theta_yz = theta_zy 
        return theta;
     }}

     // Transform Quadrupole Matrix from Cartesian to Spherical
Eigen::VectorXd PolarSite::CalculateSphericalMultipoles(const Eigen::Matrix3d& _quadrupole_cartesian){
            Eigen::Matrix3d theta = _quadrupole_cartesian ;
            Eigen::VectorXd quadrupole_polar;
            quadrupole_polar(0) = theta(2,2);
            quadrupole_polar(1) = (2./std::sqrt(3)) * theta(0,2);
            quadrupole_polar(2) = (2./std::sqrt(3)) * theta(1,2);
            quadrupole_polar(3) = (1./std::sqrt(3)) * (theta(0,0)-theta(1,1));
            quadrupole_polar(4) = (2./std::sqrt(3)) * theta(0,1);
            return quadrupole_polar;
         }

void PolarSite::Rotate(const Eigen::Matrix3d& R){
                 Eigen::Matrix3d cartesianquad = getCartesianMultipoles();
                 Eigen::Matrix3d rotated=R*cartesianquad*R.transpose();
                 CalculateSphericalMultipoles(rotated);
                 return;
             }
         

void PolarSite::Translate(const Eigen::VectorXd &shift) {
    _pos += shift;
    return;
}

void PolarSite::Induce(double wSOR) {
    // SUCCESSIVE OVERRELAXATION
    _inducedDipole_old=_inducedDipole;// Remember all previous moments
    
    _inducedDipole=(1-wSOR)*_inducedDipole_old-_Ps*(_localpermanetField+_localinducedField);
    return;  
}


/*
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

*/


}}
