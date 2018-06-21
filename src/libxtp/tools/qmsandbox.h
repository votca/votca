/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_QMSANDBOX_H
#define _VOTCA_XTP_QMSANDBOX_H

#include <stdio.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/qmtool.h>

#include <votca/xtp/qmpackagefactory.h>
#include<votca/xtp/aobasis.h>
#include<votca/xtp/aomatrix.h>

namespace votca { namespace xtp {

    
class QMSandbox : public ctp::QMTool
{
public:

    QMSandbox() { };
   ~QMSandbox() { };

    std::string Identify() { return "qmsandbox"; }

    void   Initialize(tools::Property *options);
    bool   Evaluate();


private:
    
    std::string      _orbfile;
    std::string      _mpsfiled;
    std::string      _mpsfileds;
    std::string      _mpsfileq;
    std::string      _mpsfileqs;
    
   

};

void QMSandbox::Initialize(tools::Property* options) {

    // update options with the VOTCASHARE defaults   
    //UpdateWithDefaults( options, "xtp" );
 

    std::string key = "options." + Identify();

    _orbfile = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".orbfile");
    _mpsfiled = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".dipole");
     
    _mpsfileds = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".dipole_split");
    _mpsfileq = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".quadrupole");
    _mpsfileqs = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".quadrupole_split");

    
}

bool QMSandbox::Evaluate() {
     Orbitals _orbitals;
    _orbitals.ReadFromCpt(_orbfile);
    
  #ifdef _OPENMP
 
    omp_set_num_threads(1);
    std::cout<< " Using "
                                   << omp_get_max_threads() << " threads"
                                   << std::endl;
  
#endif

    std::vector<ctp::PolarSeg*> polar_segments_dipole;
    std::vector<ctp::PolarSeg*> polar_segments_dipole_split;
    std::vector<ctp::PolarSeg*> polar_segments_quadrupole;
    std::vector<ctp::PolarSeg*> polar_segments_quadrupole_split;


    {
        vector<ctp::APolarSite*> sites = ctp::APS_FROM_MPS(_mpsfiled, 0);
        ctp::PolarSeg *newPolarSegment = new ctp::PolarSeg(0, sites);
        polar_segments_dipole.push_back(newPolarSegment);
    }   
    {
        vector<ctp::APolarSite*> sites = ctp::APS_FROM_MPS(_mpsfileds, 0);
        ctp::PolarSeg *newPolarSegment = new ctp::PolarSeg(0, sites);
        polar_segments_dipole_split.push_back(newPolarSegment);
    }        
    {
        vector<ctp::APolarSite*> sites = ctp::APS_FROM_MPS(_mpsfileq, 0);
        ctp::PolarSeg *newPolarSegment = new ctp::PolarSeg(0, sites);
        polar_segments_quadrupole.push_back(newPolarSegment);
    }        
    {
        vector<ctp::APolarSite*> sites = ctp::APS_FROM_MPS(_mpsfileqs, 0);
        ctp::PolarSeg *newPolarSegment = new ctp::PolarSeg(0, sites);
        polar_segments_quadrupole_split.push_back(newPolarSegment);
    }        
  
    Eigen::MatrixXd dmat=_orbitals.DensityMatrixGroundState();
    
    BasisSet basis;
    basis.LoadBasisSet(_orbitals.getDFTbasis());
    
    AOBasis aobasis;
    aobasis.AOBasisFill(basis,_orbitals.QMAtoms());
    
    
    AOESP esp1;
    esp1.Fillextpotential(aobasis,polar_segments_dipole_split);
    
    double edipole_split=esp1.getExternalpotential().cwiseProduct(dmat).sum();
    
    AODipole_Potential dip;
    dip.Fillextpotential(aobasis,polar_segments_dipole);
    
    double edipole=dip.getExternalpotential().cwiseProduct(dmat).sum();
    
    
    
    AOESP esp2;
    esp2.Fillextpotential(aobasis,polar_segments_quadrupole_split);
    
    double equadrupole_split=esp2.getExternalpotential().cwiseProduct(dmat).sum();
    
     AOQuadrupole_Potential quad;
    quad.Fillextpotential(aobasis,polar_segments_quadrupole);
    
    double equadrupole=quad.getExternalpotential().cwiseProduct(dmat).sum();  
    cout<<""<<endl;
    cout<<"dipole: "<<edipole<<endl;
    cout<<"dipole_split: "<<edipole_split<<endl;
    cout<<"quadrupole: "<<equadrupole<<endl;
    cout<<"quadrupole_split: "<<equadrupole_split<<endl;
    
    return true;
}






}}


#endif