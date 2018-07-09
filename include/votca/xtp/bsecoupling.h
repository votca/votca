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

#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>

#ifndef _VOTCA_XTP_BSECOUPLING_H
#define	_VOTCA_XTP_BSECOUPLING_H

namespace votca { namespace xtp {


/**
* \brief Evaluates electronic coupling elements
*
* J. Wehner,B. Baumeier, 
* JCTC DOI: 10.1021/acs.jctc.6b00935
* 
*/

class BSECoupling 
{
public:

    BSECoupling() {};
   ~BSECoupling() {};
   
   void    Initialize( tools::Property *options);
   std::string  Identify() { return "bsecoupling"; }
  
    
    Eigen::MatrixXd getJAB_singletstorage(){ return (_output_perturbation ?  JAB_singlet[0]:JAB_singlet[1]);}
       
    Eigen::MatrixXd getJAB_tripletstorage(){ return (_output_perturbation ?  JAB_triplet[0]: JAB_triplet[1]);}
    void addoutput(tools::Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB);
    
    bool CalculateCouplings(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB 
                             );  
    
  
    

     
    double getSingletCouplingElement( int levelA, int levelB, int methodindex);
    
    double getTripletCouplingElement( int levelA, int levelB, int methodindex);
   
    void setLogger( ctp::Logger* pLog ) { _pLog = pLog; }
    
private:
    
    ctp::Logger *_pLog;
  
    
    std::vector< Eigen::MatrixXd >ProjectExcitons(const Eigen::MatrixXd& _bseA_T,const Eigen::MatrixXd& _bseB_T, 
                         Eigen::MatrixXd& _H);
    
    Eigen::MatrixXd Fulldiag(const Eigen::MatrixXd& _J_dimer);
    
    Eigen::MatrixXd Perturbation(const Eigen::MatrixXd& _J_dimer);
    
    std::vector< Eigen::MatrixXd > JAB_singlet;
    std::vector< Eigen::MatrixXd > JAB_triplet;

    bool _doTriplets;
    bool _doSinglets;
    bool _output_perturbation;
    int _levA;
    int _levB;
    int _occA;
    int _unoccA;
    int _occB;
    int _unoccB;
    double      _degeneracy;
    int         _openmp_threads;
    
    
    //_levA und _bseA_exc are the same but the stupid int unsigned conversion stuff, so leave it for now
    
    unsigned _bse_exc;
    
    unsigned _bseA_exc;
    unsigned _bseB_exc;
    
    unsigned _ct;
    
     Eigen::MatrixXd ctAB;
     Eigen::MatrixXd ctBA;
     Eigen::MatrixXd _kap;
     Eigen::MatrixXd _kbp;
    
    
};

}}

#endif	/* _VOTCA_XTP_BSECOUPLING_H */


