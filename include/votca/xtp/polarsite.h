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

#include "qmatom.h"
namespace votca { namespace xtp {
    /**
    \brief Class to represent Atom/Site in electrostatic+polarisation 

     The units are atomic units, e.g. Bohr, Hartree.By default a PolarSite cannot be polarised.
*/
class PolarSite
{

public:

    PolarSite(int id, const std::string& element, const Eigen::Vector3d& pos);
            
    PolarSite(int id, const std::string& element)
                    :PolarSite(id,element,Eigen::Vector3d::Zero()){
                };

    PolarSite(const QMAtom& atom):PolarSite(atom.getAtomID(),atom.getElement(),atom.getPos()){
        Eigen::VectorXd m=Eigen::VectorXd::Zero(1);
        m<<atom.getNuccharge();
        setMultipole(m);
    }
      

    int getId() const{ return _id; }
    int getRank()const{return _rank;}
    const std::string &getElement() const{ return _element; }
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
    
    // COORDINATES TRANSFORMATION
    void Translate(const Eigen::VectorXd &shift);
    void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos);
 
    // MULTIPOLES DEFINITION
    
    double getCharge() const{return _multipole(0);}
    const Eigen::VectorXd& getPermMultipole()const {return _multipole;}//Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s,...[NOT following Stone order for dipoles]
    
    Eigen::Vector3d getDipole()const{
        Eigen::Vector3d dipole=Eigen::Vector3d::Zero();
        if(_isPolarisable){
            dipole+=_inducedDipole;
        }
        if(_rank>0){
            dipole+=_multipole.segment<3>(1);
        }
        return dipole;
    }
    
    Eigen::Matrix3d CalculateCartesianMultipole()const; 
    static Eigen::VectorXd CalculateSphericalMultipole(const Eigen::Matrix3d& _quadrupole_cartesian);
    
    Eigen::Vector3d getField()const{return _localpermanetField+_localinducedField;}
    
    double getPotential()const{return PhiP+PhiU;}
    
    void WriteMpsLine(std::ostream &out, const std::string& unit = "bohr")const;
    void Induce(double wSOR);
       
    double InteractStatic(PolarSite& otherSite);
    
    double InteractInduction(PolarSite& otherSite, double a=0.39);
    
    double InductionWork() const{ return -0.5*_inducedDipole.transpose()*getField();}
    
    void WriteToCpt(CptLoc parent)const;

   void ReadFromCpt(CptLoc parent);
    
    
private:
       
    void calcRank(); 
    Eigen::MatrixXd FillTholeInteraction(const PolarSite& otherSite, double a);
    Eigen::MatrixXd FillInteraction(const PolarSite& otherSite);
    
    int     _id;
    std::string  _element;
    Eigen::Vector3d _pos;
    int     _rank;
    Eigen::VectorXd _multipole; //Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s
    
    //required for polarisation
    bool _isPolarisable;
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
