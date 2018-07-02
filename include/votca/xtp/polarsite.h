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

#ifndef __VOTCA_XTP_POLARSITE_H
#define __VOTCA_XTP_POLARSITE_H

#include <votca/xtp/eigen.h>
#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/segment.h>
namespace votca { namespace xtp {
    
class PolarSite
{

public:
    
    PolarSite(int id, const std::string name)
            : _id(id), _name(name),
            PhiP(0.0),PhiU(0.0){_pos=Eigen::Vector3d::Zero(3);}

    PolarSite(int id, const std::string name, const Eigen::Vector3d pos)
            : _id(id), _name(name),
            PhiP(0.0),PhiU(0.0){_pos=pos;}
        
     
    // GET & SET & IMPORT FUNCTIONS
    int    getId() const{ return _id; }
    int getRank()const{return _rank;}
    const std::string &getName() const{ return _name; }
    const Eigen::Vector3d &getPos() const{ return _pos; }
    
    
    void setMultipoles(const Eigen::VectorXd& multipoles){
        _multipoles=multipoles;
        calcRank();
    }
   
    void setSegment(Segment* seg){
        _segment=seg;
    }
    void setFragment(Fragment* frag){
        _fragment=frag;
    }
    
    void setTopology(Topology* top){
        _top=top;
    }
    // COORDINATES TRANSFORMATION
    void    Translate(const Eigen::VectorXd &shift);
    void Rotate(const Eigen::Matrix3d& R);
 
    // MULTIPOLES DEFINITION
    const Eigen::VectorXd & getMultipoles()const {return _multipoles;}//Q00,Q10, Q11s, Q11c,Q20, Q21c, Q21s, Q22c, Q22s,...[following Stone order]
    //Eigen::MatrixXd getQuad_cartesian(){return _quadrupole_cartesian;} //theta matrix
    Eigen::VectorXd getQuad()const{return _quadrupole_polar;}//Q20, Q21c, Q21s, Q22c, Q22s
    Eigen::Vector3d getCartesianDipoles() const;// mu_x,mu_y,mu_z  
    Eigen::Matrix3d getCartesianMultipoles(); 
    Eigen::VectorXd CalculateSphericalMultipoles(const Eigen::Matrix3d& _quadrupole_cartesian);
    

    //void WriteMpsLine(std::ostream &out, string unit = "angstrom");
    void Induce(double wSOR);
    
    Eigen::MatrixXd Interaction(const PolarSite& otherSite);
    
    
private:
 void calcRank(); 
    int     _id;
    std::string  _name;
    Eigen::Vector3d _pos;
    int     _rank;

    Eigen::Matrix3d _Ps;
    Eigen::Vector3d _localpermanetField;
    Eigen::Vector3d _localinducedField;
    Eigen::Vector3d _inducedDipole;
    Eigen::Vector3d _inducedDipole_old;
    
    Eigen::VectorXd _multipoles; //Q00,Q10, Q11c, Q11s,Q20, Q21c, Q21s, Q22c, Q22s
    Eigen::Vector3d _dipole; // mu_x, mu_y, mu_z
    Eigen::Matrix3d _quadrupole_cartesian; //theta matrix
    Eigen::VectorXd _quadrupole_polar; // Q20, Q21c, Q21s, Q22c, Q22s
    
    Segment* _segment;
    Fragment* _fragment;
    Topology* _top;
    double PhiP;                            // Electric potential (due to perm.)
    double PhiU;                            // Electric potential (due to indu.)

};


}}


#endif
