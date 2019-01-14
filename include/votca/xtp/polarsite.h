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
#include <votca/xtp/staticsite.h>

namespace votca { namespace xtp {\
    
  
    /**
    \brief Class to represent Atom/Site in electrostatic+polarisation 

     The units are atomic units, e.g. Bohr, Hartree.By default a PolarSite cannot be polarised.
*/
class PolarSite : public StaticSite
{

public:
    

    PolarSite(int id, std::string element, Eigen::Vector3d pos);
            
    PolarSite(int id, std::string element)
                    :PolarSite(id,element,Eigen::Vector3d::Zero()){
                };

    PolarSite(const CheckpointReader& r){
        ReadFromCpt(r);
    }

    PolarSite(const QMAtom& atom, double charge):PolarSite(atom.getAtomID(),atom.getElement(),atom.getPos()){
        setCharge(charge);
    }
      
        
    void setPolarisation(const Eigen::Matrix3d pol){
        _Ps=pol;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
        es.computeDirect(_Ps,Eigen::EigenvaluesOnly);
        _eigendamp=es.eigenvalues().maxCoeff();
    }
    
    void ResetInduction(){
        PhiU=0.0;
        _inducedDipole=Eigen::Vector3d::Zero();
        _inducedDipole_old=Eigen::Vector3d::Zero();
        _localinducedField=Eigen::Vector3d::Zero();
    }
     
    // MULTIPOLES DEFINITION
       
    Eigen::Vector3d getDipole()const{
        Eigen::Vector3d dipole=Eigen::Vector3d::Zero();
        dipole+=_inducedDipole;
        dipole+=_multipole.segment<3>(1);
        return dipole;
    }
   
    void Induce(double wSOR);
    
    double InductionWork() const{ return -0.5*_inducedDipole.transpose()*getField();}
    
    static std::string Identify(){return "polarsite";}
    
private:
    
    int     _id;
    std::string  _element;
    Eigen::Vector3d _pos;
    int     _rank;

    Vector9d _multipole; //Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s
    
    
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
