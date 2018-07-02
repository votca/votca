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
    //It calculates the rank in the spherical case
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
Eigen::Vector3d PolarSite::getCartesianDipoles()const{
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

Eigen::MatrixXd PolarSite::Interaction(const PolarSite& otherSite){
    
    const Eigen::Vector3d& posB= otherSite.getPos();
    const Eigen::Vector3d& posA=this->getPos();
    
    
    Eigen::MatrixXd interaction=Eigen::MatrixXd::Zero(this->getMultipoles().size(),otherSite.getMultipoles().size());
    Eigen::Vector3d r_AB = posB-posA; //Vector of the distance between polar sites
    Eigen::Vector3d e_AB=r_AB.normalized(); //Unit Vector that points from a site to b site
    Eigen::Vector3d pos_a=e_AB; //unit vector on the sites reciprocal direction; This points toward A
    Eigen::Vector3d pos_b=-e_AB; //unit vector on the sites reciprocal direction; This points toward B   
    double R = r_AB.norm(); //Norm of distance vector
    double fac0=1/R;
    Eigen::Matrix3d c = Eigen::Matrix3d::Identity(); //Following Stone's notation. Here we are assuming that the sites are in the same localframe (Jens dixit)
    
    
    int rank = _rank;
    
    int rankA=this->getRank();
    int rankB=otherSite.getRank();
    
    double fac1=std::pow(fac0,2);
    double fac2=std::pow(fac0,3);
    double fac3=std::pow(fac0,4);
    double fac4=std::pow(fac0,5);
    
    //Charge-Charge Interaction
    interaction(0,0)=fac0;
    if ( rankA>0 ){
        //Dipole-Charge Interaction
        interaction.block(3,1,1,0)=fac1*pos_a;
        if ( rankB>0 ){
           //Charge-Dipole Interaction
            interaction.block(1,3,0,1)=fac1*pos_b;
            //Dipole-Dipole Interaction 
            interaction.block(3,3,1,1)=fac2*(3*pos_a*pos_b.transpose())+c; 
            }
    }
    if ( rankA > 1 ){
        interaction(0,4)=fac2*0.5*(3*std::pow(pos_b(2),2)-1); //T00,20
        interaction(0,5)=fac2*std::sqrt(3)*pos_b(0)*pos_b(2); //T00,21c
        interaction(0,6)=fac2*std::sqrt(3)*pos_b(1)*pos_b(2); //T00,21s
        interaction(0,7)=fac2*0.5*std::sqrt(3)*(std::pow(pos_b(0),2)-std::pow(pos_b(1),2)); //T00,22c
        interaction(0,8)=fac2*std::sqrt(3)*pos_b(0)*pos_b(1); //T00,22s
        
        interaction.block(3,1,1,4)=0.5*fac3*(15*std::pow(pos_b(2),2)*pos_a+6*pos_b(2)*c.block(3,1,0,2)-3*pos_a); //dipole(order mu_x,mu_y,mu_z)-20
        interaction.block(3,1,1,5)=fac3*std::sqrt(3)*(pos_b(1)*c.block(3,1,0,2)+c.block(3,1,0,0)*pos_b(2)+5*pos_b(0)*pos_b(2)*pos_a);//dipole(order mu_x,mu_y,mu_z)-21c
        interaction.block(3,1,1,6)=fac3*std::sqrt(3)*(pos_b(1)*c.block(3,1,0,2)+c.block(3,1,0,1)*pos_b(2)+5*pos_b(1)*pos_b(2)*pos_a);//dipole(order mu_x,mu_y,mu_z)-21s
        interaction.block(3,1,1,7)=fac3*0.5*std::sqrt(3)*(5*(std::pow(pos_b(0),2)-std::pow(pos_b(1),2))*pos_a+2*pos_b(0)*c.block(3,1,0,0)-2*pos_b(1)*c.block(3,1,0,1));//dipole(order mu_x,mu_y,mu_z)-22c
        interaction.block(3,1,1,8)=fac3*std::sqrt(3)*(5*pos_b(0)*pos_b(1)*pos_a+pos_b(0)*c.block(1,3,1,0)+pos_b(1)*c.block(3,1,0,0));//dipole(order mu_x,mu_y,mu_z)-22s
        if (rankB > 1 ){
            interaction(4,0)=fac2*0.5*(3*std::pow(pos_a(2),2)-1); //T20,00
            interaction(5,0)=fac2*std::sqrt(3)*pos_a(0)*pos_a(2); //21c,00
            interaction(6,0)=fac2*std::sqrt(3)*pos_a(1)*pos_a(2); //T21s,000
            interaction(7,0)=fac2*0.5*std::sqrt(3)*(std::pow(pos_a(0),2)-std::pow(pos_a(1),2)); //T22c,00
            interaction(8,0)=fac2*std::sqrt(3)*pos_a(0)*pos_a(1); //T00,22s
            interaction.block(3,1,1,4)=0.5*fac3*(15*std::pow(pos_b(2),2)*pos_a+6*pos_b(2)*c.block(3,1,0,2)-3*pos_a); //dipole(order mu_x,mu_y,mu_z)-20
            interaction.block(3,1,4,1)=fac3*std::sqrt(3)*(pos_b(1)*c.block(3,1,0,2)+c.block(3,1,0,0)*pos_b(2)+5*pos_b(0)*pos_b(2)*pos_a);//dipole(order mu_x,mu_y,mu_z)-21c
            interaction.block(3,1,1,6)=fac3*std::sqrt(3)*(pos_b(1)*c.block(3,1,0,2)+c.block(3,1,0,1)*pos_b(2)+5*pos_b(1)*pos_b(2)*pos_a);//dipole(order mu_x,mu_y,mu_z)-21s
            interaction.block(3,1,1,7)=fac3*0.5*std::sqrt(3)*(5*(std::pow(pos_b(0),2)-std::pow(pos_b(1),2))*pos_a+2*pos_b(0)*c.block(3,1,0,0)-2*pos_b(1)*c.block(3,1,0,1));//dipole(order mu_x,mu_y,mu_z)-22c
            interaction.block(3,1,8,1)=fac3*std::sqrt(3)*(5*pos_b(0)*pos_b(1)*pos_a+pos_b(0)*c.block(1,3,1,0)+pos_b(1)*c.block(3,1,0,0));//dipole(order mu_x,mu_y,mu_z)-22s
            
            interaction(4,4)=fac4*(3./4.)*(35*std::pow(pos_a(2)*pos_b(2),2)-5*std::pow(pos_a(2),2)
                -5*std::pow(pos_b(2),2)+20*pos_a(2)*pos_b(2)*c(2,2)+2*std::pow(c(2,2),2)+1); //T20,20
            interaction(4,5)=0.5*fac4*std::sqrt(3)*(35*std::pow(pos_a(2),2)*pos_b(0)*pos_b(2)
                -5*pos_b(0)*pos_b(2)+10*pos_a(2)*pos_b(0)*c(2,2)+10*pos_a(2)*pos_b(2)*c(2,1)+2*c(2,0)*c(2,2)); //T20,21c
            interaction(4,6)=0.5*fac4*std::sqrt(3)*(35*std::pow(pos_a(2),2)*pos_b(1)*pos_b(2)
                -5*pos_b(1)*pos_b(2)+10*pos_a(2)*pos_b(1)*c(2,2)+10*pos_a(2)*pos_b(2)*c(2,1)+2*c(2,1)*c(2,2)) ; //T20,21s
            interaction(4,7)=0.25*fac4*std::sqrt(3)*(35*std::pow(pos_a(2)*pos_b(0),2)-35*std::pow(pos_a(2)*pos_b(1),2)-5*std::pow(pos_b(0),2)
                +5*std::pow(pos_b(1),2)+20*pos_a(2)*pos_b(0)*c(2,0)-20*pos_a(2)*pos_b(1)*c(2,1)+2*std::pow(c(2,0),2)-2*std::pow(c(2,1),2)); //T20,22c
            interaction(4,8)=0.5*fac4*std::sqrt(3)*(35*std::pow(pos_a(2),2)*pos_b(0)*pos_b(1)
                -5*pos_b(0)*pos_b(1)+10*pos_a(2)*pos_b(0)*c(2,1)+10*pos_a(2)*pos_b(1)*c(2,0)+2*c(2,0)*c(2,1)) ; //T20,22s
            interaction(5,5)=fac4*(35*pos_a(0)*pos_a(2)*pos_b(0)*pos_b(2)+5*pos_a(0)*pos_b(0)*c(2,2)+5*pos_a(0)*pos_b(2)*c(2,0)+
                5*pos_a(2)*pos_b(0)*c(0,2)+5*pos_a(2)*pos_b(2)*c(0,0)+c(0,0)*c(2,2)+c(0,2)*c(2,0)) ; //T21c,21c
            interaction(5,6)=fac4*(35*pos_a(0)*pos_a(2)*pos_b(1)*pos_b(2)+5*pos_a(0)*pos_b(1)*c(2,2)+5*pos_a(0)*pos_b(2)*c(2,1)+
                5*pos_a(2)*pos_b(1)*c(0,2)+5*pos_a(2)*pos_b(2)*c(0,1)+c(0,1)*c(2,2)+c(0,2)*c(2,1)) ; //T21c,21s
            interaction(5,7)=0.5*fac4*(35*pos_a(0)*pos_a(2)*std::pow(pos_b(0),2)-35*pos_a(0)*pos_a(2)*std::pow(pos_b(1),2)
                +10*pos_a(0)*pos_b(0)*c(2,0)-10*pos_a(0)*pos_b(1)*c(2,1)+10*pos_a(0)*pos_b(0)*c(0,0)-10*pos_a(2)*pos_b(1)*c(0,1)+2*c(0,0)*c(2,0)-2*c(0,1)*c(2,1)); //T21c,22c
            interaction(5,8)=fac4*(35*pos_a(0)*pos_a(2)*pos_b(0)*pos_b(1)+5*pos_a(0)*pos_b(0)*c(2,1)+5*pos_a(0)*pos_b(1)*c(2,0)+
                5*pos_a(2)*pos_b(0)*c(0,1)+5*pos_a(2)*pos_b(1)*c(0,0)+c(0,0)*c(2,1)+c(0,1)*c(2,0)) ; //T21c,22s
            interaction(6,6)=fac4*(35*pos_a(1)*pos_a(2)*pos_b(1)*pos_b(2)+5*pos_a(1)*pos_b(1)*c(2,2)+5*pos_a(1)*pos_b(2)*c(2,1)+
                5*pos_a(2)*pos_b(1)*c(1,2)+5*pos_a(2)*pos_b(2)*c(1,1)+c(1,1)*c(2,2)+c(1,2)*c(2,1)) ; //T21s,21s
            interaction(6,7)=0.5*fac4*(35*pos_a(1)*pos_a(2)*std::pow(pos_b(0),2)-35*pos_a(1)*pos_a(2)*std::pow(pos_b(1),2)
                +10*pos_a(1)*pos_b(2)*c(2,0)-10*pos_a(1)*pos_b(1)*c(2,1)+10*pos_a(2)*pos_b(0)*c(1,0)-10*pos_a(2)*pos_b(1)*c(1,1)+2*c(1,0)*c(2,0)-2*c(1,1)*c(2,1)); //T21s,22c
            interaction(6,8)=fac4*(35*pos_a(1)*pos_a(2)*pos_b(0)*pos_b(1)+5*pos_a(1)*pos_b(0)*c(2,1)+5*pos_a(1)*pos_b(1)*c(2,1)+
                5*pos_a(2)*pos_b(0)*c(1,1)+5*pos_a(2)*pos_b(1)*c(1,2)+c(1,0)*c(2,1)+c(1,1)*c(2,0)) ; //T21s,22s
            interaction(7,7)=0.25*fac4*(35*std::pow(pos_a(0)*pos_b(0),2)-35*std::pow(pos_a(0)*pos_b(1),2)-35*std::pow(pos_a(1)*pos_b(0),2)
                +35*std::pow(pos_a(1)*pos_b(1),2)+20*pos_a(0)*pos_b(0)*c(0,0)-20*pos_a(0)*pos_b(1)*c(0,1)-20*pos_a(1)*pos_b(0)*c(1,0)
                +20*pos_a(0)*pos_b(0)*c(1,1)+2*std::pow(c(0,0),2)-2*std::pow(c(0,1),2)-2*std::pow(c(1,0),2)+2*std::pow(c(1,1),2)); //T22c,22c
            interaction(7,8)= 0.5*fac4*(35*pos_b(0)*pos_b(1)*std::pow(pos_a(0),2)-35*pos_b(1)*pos_b(2)*std::pow(pos_a(1),2)
                +10*pos_a(0)*pos_b(0)*c(0,1)+10*pos_a(0)*pos_b(1)*c(0,0)-10*pos_a(1)*pos_b(0)*c(1,1)-10*pos_a(1)*pos_b(1)*c(1,2)+2*c(0,0)*c(0,1)-2*c(1,0)*c(1,1)); //T22c,22s
            interaction(8,8)=0.5*fac4*(35*pos_a(0)*pos_a(1)*pos_b(0)*pos_b(1)+5*pos_a(0)*pos_b(0)*c(1,1)+5*pos_a(0)*pos_b(1)*c(1,0)+5*pos_a(1)*pos_b(0)*c(0,1)
                +5*pos_a(1)*pos_b(1)*c(0,0)+c(0,0)*c(1,1)+c(0,1)*c(1,0)); //T22s,22s
            interaction(5,4)=0.5*fac4*std::sqrt(3)*(35*std::pow(pos_b(2),2)*pos_a(0)*pos_a(2)
                -5*pos_a(0)*pos_a(2)+10*pos_b(2)*pos_a(0)*c(2,2)+10*pos_b(2)*pos_a(2)*c(1,2)+2*c(0,2)*c(2,2)); //T21c,20
            interaction(6,4)=0.5*fac4*std::sqrt(3)*(35*std::pow(pos_b(2),2)*pos_a(1)*pos_a(2)
                -5*pos_a(1)*pos_a(2)+10*pos_b(2)*pos_a(1)*c(2,2)+10*pos_b(2)*pos_a(2)*c(1,2)+2*c(1,2)*c(2,2)) ; //T21s,20
            interaction(6,5)=fac4*(35*pos_b(0)*pos_b(2)*pos_a(1)*pos_a(2)+5*pos_b(0)*pos_a(1)*c(2,2)+5*pos_b(0)*pos_a(2)*c(1,2)+
                5*pos_b(2)*pos_a(1)*c(2,0)+5*pos_b(2)*pos_a(2)*c(1,0)+c(1,0)*c(2,2)+c(2,0)*c(1,2)) ; //T21s,21c
            interaction(7,4)=0.25*fac4*std::sqrt(3)*(35*std::pow(pos_b(2)*pos_a(0),2)-35*std::pow(pos_b(2)*pos_a(1),2)-5*std::pow(pos_a(0),2)
                +5*std::pow(pos_a(1),2)+20*pos_b(2)*pos_a(0)*c(0,2)-20*pos_b(2)*pos_a(1)*c(1,2)+2*std::pow(c(0,2),2)-2*std::pow(c(1,2),2)); //T22c,20
            interaction(7,5)=0.5*fac4*(35*pos_b(0)*pos_b(2)*std::pow(pos_a(0),2)-35*pos_b(0)*pos_b(2)*std::pow(pos_a(1),2)
                +10*pos_b(0)*pos_a(0)*c(0,2)-10*pos_b(0)*pos_a(1)*c(1,2)+10*pos_b(0)*pos_a(0)*c(0,0)-10*pos_b(2)*pos_a(1)*c(1,0)+2*c(0,0)*c(0,2)-2*c(1,0)*c(1,2)); //T22c,21c
            interaction(7,6)=0.5*fac4*(35*pos_b(1)*pos_b(2)*std::pow(pos_a(0),2)-35*pos_b(1)*pos_b(2)*std::pow(pos_a(1),2)
                +10*pos_b(1)*pos_a(2)*c(0,2)-10*pos_b(1)*pos_a(1)*c(1,2)+10*pos_b(2)*pos_a(0)*c(0,1)-10*pos_b(2)*pos_a(1)*c(1,1)+2*c(0,1)*c(0,2)-2*c(1,1)*c(1,2)); //T22c,21s
            interaction(8,4)=0.5*fac4*std::sqrt(3)*(35*std::pow(pos_b(2),2)*pos_a(0)*pos_a(1)
                -5*pos_a(0)*pos_a(1)+10*pos_b(2)*pos_a(0)*c(1,2)+10*pos_b(2)*pos_a(1)*c(0,2)+2*c(0,2)*c(1,2)) ; //T22s,20
            interaction(8,5)=fac4*(35*pos_b(0)*pos_b(2)*pos_a(0)*pos_a(1)+5*pos_b(0)*pos_a(0)*c(1,2)+5*pos_b(0)*pos_a(1)*c(0,2)+
                5*pos_b(2)*pos_a(0)*c(1,0)+5*pos_b(2)*pos_a(1)*c(0,0)+c(0,0)*c(1,2)+c(1,0)*c(0,2)) ; //T22s,21c
            interaction(8,6)=fac4*(35*pos_b(1)*pos_b(2)*pos_a(0)*pos_a(1)+5*pos_b(1)*pos_a(0)*c(1,2)+5*pos_b(1)*pos_a(1)*c(1,2)+
                5*pos_b(2)*pos_a(0)*c(1,1)+5*pos_b(2)*pos_a(1)*c(2,1)+c(0,1)*c(1,2)+c(1,1)*c(0,2)) ; //T22s,21s
            interaction(8,7)= 0.5*fac4*(35*pos_a(0)*pos_a(1)*std::pow(pos_b(0),2)-35*pos_a(1)*pos_a(2)*std::pow(pos_b(1),2)
                +10*pos_b(0)*pos_a(0)*c(1,0)+10*pos_b(0)*pos_a(1)*c(0,0)-10*pos_b(1)*pos_a(0)*c(1,1)-10*pos_b(1)*pos_a(1)*c(2,1)+2*c(0,0)*c(1,0)-2*c(0,1)*c(1,1)); //T22s,22c
             
        }
    }
    
        return interaction; //in units of 4piepsilon0
        
}
    
    
       





/*
void PolarSite::CalculateMultipolesEnergy() {
    Eigen::VectorXd mp_a; //Multipoles on site a
    Eigen::VectorXd mp_b; //Multipoles on site b
    
}
 */

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
