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

#ifndef __VOTCA_XTP_STATICSITE_H
#define __VOTCA_XTP_STATICSITE_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/qmatom.h>

namespace votca { namespace xtp {\
    
  
    /**
    \brief Class to represent Atom/Site in electrostatic 

     The units are atomic units, e.g. Bohr, Hartree.
*/
class StaticSite
{

public:
    

    StaticSite(int id, std::string element, Eigen::Vector3d pos);
            
    StaticSite(int id, std::string element)
                    :StaticSite(id,element,Eigen::Vector3d::Zero()){
                };

    StaticSite(const CheckpointReader& r){
        ReadFromCpt(r);
    }

    StaticSite(const QMAtom& atom, double charge):StaticSite(atom.getAtomID(),atom.getElement(),atom.getPos()){
        setCharge(charge);
    }
      

    int getId() const{ return _id; }
    int getRank()const{return _rank;}
    const std::string &getElement() const{ return _element; }
    const Eigen::Vector3d &getPos() const{ return _pos; }
       
    
    void setMultipole(const Vector9d& multipole){
        _multipole=multipole;
        calcRank();
    }

    void setCharge(double q){
        _multipole(0)=q;
        calcRank();
    }
    
    
    // COORDINATES TRANSFORMATION
    void Translate(const Eigen::VectorXd &shift);
    void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos);
 
    // MULTIPOLES DEFINITION
    
    double getCharge() const{return _multipole(0);}
    const Vector9d& getPermMultipole()const {return _multipole;}//Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s,...[NOT following Stone order for dipoles]
    
    
    virtual Eigen::Vector3d getDipole()const{
        return _multipole.segment<3>(1);
    }
   
    Eigen::Matrix3d CalculateCartesianMultipole()const; 
    static Eigen::VectorXd CalculateSphericalMultipole(const Eigen::Matrix3d& quadrupole_cartesian);
    
    Eigen::Vector3d getField()const{return _localpermanetField;}
    
    double getPotential()const{return PhiU;}
    
    std::string WriteMpsLine(std::string unit = "bohr")const;
    
    void WriteToCpt(const CheckpointWriter& w)const;

    void ReadFromCpt(const CheckpointReader& r);
   
    virtual std::string Identify(){return "staticsite";}
    
    
protected:
       
    void calcRank(); 
    int     _id;
    std::string  _element;
    Eigen::Vector3d _pos;
    int     _rank;

    Eigen::VectorXd _multipole; //Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s

    
    Eigen::Matrix3d _Ps;
    Eigen::Vector3d _localpermanetField;
    Eigen::Vector3d _localinducedField;
    Eigen::Vector3d _inducedDipole;
    Eigen::Vector3d _inducedDipole_old;
    double _eigendamp;
    
    double PhiP;                            // Electric potential (due to perm.)
    double PhiU;                            // Electric potential (due to indu.)
    
};


}}


#endif
