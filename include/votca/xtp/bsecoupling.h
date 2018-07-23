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

#include <votca/xtp/couplingbase.h>

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

class BSECoupling : CouplingBase
{
public:
   
   void    Initialize( tools::Property &options);
   std::string  Identify() { return "bsecoupling"; }
  
    
    Eigen::MatrixXd getJAB_singletstorage(){ return (_output_perturbation ?  JAB_singlet[0]:JAB_singlet[1]);}
       
    Eigen::MatrixXd getJAB_tripletstorage(){ return (_output_perturbation ?  JAB_triplet[0]: JAB_triplet[1]);}
    void Addoutput(tools::Property & type_summary,const Orbitals& orbitalsA, 
                               const Orbitals& orbitalsB);
    

    
    void CalculateCouplings(   const Orbitals& orbitalsA, 
                               const Orbitals& orbitalsB, 
                               const Orbitals& orbitalsAB 
                             );
     
   
    
    
private:
    
    void WriteToProperty(const Orbitals& orbitalsA, const Orbitals& orbitalsB,
                        tools::Property& summary, int stateA, int stateB, double JAB);
    
    double getSingletCouplingElement( int levelA, int levelB, int methodindex);
    
    double getTripletCouplingElement( int levelA, int levelB, int methodindex);
    
   
  
    
    std::vector< Eigen::MatrixXd >ProjectExcitons(const Eigen::MatrixXd& bseA_T,const Eigen::MatrixXd& bseB_T, 
                         Eigen::MatrixXd& H);
    
    Eigen::MatrixXd Fulldiag(const Eigen::MatrixXd& J_dimer);
    
    Eigen::MatrixXd Perturbation(const Eigen::MatrixXd& J_dimer);
    
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
    int         _openmp_threads;
    
    
    int _bse_exc;
    
    int _ct;
    
     Eigen::MatrixXd ctAB;
     Eigen::MatrixXd ctBA;
     Eigen::MatrixXd _kap;
     Eigen::MatrixXd _kbp;
    
    
};

}}

#endif	/* _VOTCA_XTP_BSECOUPLING_H */


