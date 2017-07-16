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

#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>

#ifndef _VOTCA_XTP_BSECOUPLING_H
#define	_VOTCA_XTP_BSECOUPLING_H

namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;
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
   
   void    Initialize( Property *options);
    string  Identify() { return "bsecoupling"; }
    bool get_doSinglets(){ return _doSinglets;}
    bool get_doTriplets(){ return _doTriplets;}
    
    ub::matrix<double> getJAB_singletstorage(){return JAB_singlet;}
    ub::matrix<double> getJAB_tripletstorage(){return JAB_triplet;}
    void addoutput(Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB);
    
    bool CalculateCouplings(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB 
                             );  
    
  
    

     
    double getSingletCouplingElement( int levelA, int levelB);
    
    double getTripletCouplingElement( int levelA, int levelB);
    double getSingletDimerEnergy( int level);
    double getTripletDimerEnergy( int level);
    void setLogger( ctp::Logger* pLog ) { _pLog = pLog; }
    
private:
    
    ctp::Logger *_pLog;
  
    
    bool ProjectExcitons(const ub::matrix<double>& _bseA_T,const ub::matrix<double>& _bseB_T, 
                         ub::matrix<double>& _H, ub::matrix<double>& _J );
    
    ub::matrix<double> JAB_singlet;
    ub::matrix<double> JAB_triplet;

    bool _doTriplets;
    bool _doSinglets;
    bool _do_perturbation;
    bool _do_full_diag;
    int _levA;
    int _levB;
    int _occA;
    int _unoccA;
    int _occB;
    int _unoccB;
    double      _degeneracy;
    int         _openmp_threads;
    
    
     ub::matrix<double> ctAB;
     ub::matrix<double> ctBA;
     ub::matrix<double> _kap;
     ub::matrix<double> _kbp;
    
    
};

}}

#endif	/* _VOTCA_XTP_BSECOUPLING_H */


