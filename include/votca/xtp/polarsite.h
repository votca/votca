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
#include <votca/xtp/topology.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/segment.h>
namespace votca { namespace xtp {
    /**
    \brief Class to represent Atom/Site in electrostatic+polarisation 

     The units are atomic units, e.g. Bohr, Hartree.By default a PolarSite cannot be polarised.
*/
class PolarSite
{

public:
    


    PolarSite(int id, const std::string& name, const Eigen::Vector3d& pos)
            : _id(id), _name(name),_isPolarisable(false),
            PhiP(0.0),PhiU(0.0),
            _quadrupole(_multipole,4,5),_dipole(_multipole,1,3),
            _pos(pos),
            _localpermanetField(Eigen::Vector3d::Zero()),
            _localinducedField(Eigen::Vector3d::Zero()){};
            
    PolarSite(int id, const std::string& name)
         : _id(id), _name(name),_isPolarisable(false),
         PhiP(0.0),PhiU(0.0),
         _quadrupole(_multipole,4,5),
         _dipole(_multipole,1,3),
         _pos(Eigen::Vector3d::Zero()),
        _localpermanetField(Eigen::Vector3d::Zero()),
        _localinducedField(Eigen::Vector3d::Zero()){};
        
    
        
     
    int getId() const{ return _id; }
    int getRank()const{return _rank;}
    const std::string &getName() const{ return _name; }
    const Eigen::Vector3d &getPos() const{ return _pos; }
    
    bool isPolarisable() const{ return _isPolarisable;}
    
    void setPolarisable(bool polarisable){
        _isPolarisable=polarisable;
    }
    
    void setMultipole(const Eigen::VectorXd& multipole){
        _multipole=multipole;
        calcRank();
    }
    
    void setPolarisation(const Eigen::Matrix3d pol){
        _Ps=pol;
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
    void Translate(const Eigen::VectorXd &shift);
    void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos);
 
    // MULTIPOLES DEFINITION
    
    double getCharge() const{return _multipole(0);}
    const Eigen::VectorXd& getPermMultipole()const {return _multipole;}//Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s,...[NOT following Stone order for dipoles]
    const Eigen::VectorBlock< const Eigen::VectorXd,5 >& getQuadrupole()const{return _quadrupole;}//Q20, Q21c, Q21s, Q22c, Q22s
    const Eigen::VectorBlock< const Eigen::VectorXd,3 >& getPermDipole() const{return _dipole;}// mu_x,mu_y,mu_z  
    const Eigen::Vector3d& getInducedDipole()const{return _inducedDipole;}
    Eigen::Matrix3d CalculateCartesianMultipole(); 
    static Eigen::VectorXd CalculateSphericalMultipole(const Eigen::Matrix3d& _quadrupole_cartesian);
    
    Eigen::Vector3d getField(){return _localpermanetField+_localinducedField;}
    
    double getPotential(){return PhiP+PhiU;}
    
    void WriteMpsLine(std::ostream &out, const std::string& unit = "bohr");
    void Induce(double wSOR);
       
    double InteractStatic(PolarSite& otherSite);
    
    double InteractInduction(PolarSite& otherSite);
    
private:
    
    
    void calcRank(); 
    Eigen::MatrixXd FillInteraction(const PolarSite& otherSite);
    
    int     _id;
    std::string  _name;
    Eigen::Vector3d _pos;
    int     _rank;

    
    bool _isPolarisable;
    Eigen::Matrix3d _Ps;
    Eigen::Vector3d _localpermanetField;
    Eigen::Vector3d _localinducedField;
    Eigen::Vector3d _inducedDipole;
    Eigen::Vector3d _inducedDipole_old;
    
    Eigen::VectorXd _multipole; //Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s
    Eigen::VectorBlock< const Eigen::VectorXd,5 > _quadrupole;
    Eigen::VectorBlock< const Eigen::VectorXd,3 > _dipole;
    
    
    
    
    Segment* _segment;
    Fragment* _fragment;
    Topology* _top;
    double PhiP;                            // Electric potential (due to perm.)
    double PhiU;                            // Electric potential (due to indu.)

};


}}


#endif
