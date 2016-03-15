/*
 *            Copyright 2009-2016 The VOTCA Development Team
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
#include <votca/xtp/logger.h>

#ifndef _VOTCA_XTP_BSECOUPLING_H
#define	_VOTCA_XTP_BSECOUPLING_H

namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;
/**
* \brief Evaluates electronic coupling elements
*
* B. Baumeier, J. Kirkpatrick, D. Andrienko, 
* Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
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
    
    ub::matrix<float> getJAB_singletstorage(){return JAB_singlet;}
    ub::matrix<float> getJAB_tripletstorage(){return JAB_triplet;}
    void addoutput(Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB);
    
    bool CalculateCouplings(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB 
                             );  
    
    bool CalculateCouplings_OLD(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB, 
                               ub::matrix<float>* _JAB_singlet);  
    
    

     
    float getSingletCouplingElement( int levelA, int levelB);
    
    float getTripletCouplingElement( int levelA, int levelB);
    float getSingletDimerEnergy( int level);
    float getTripletDimerEnergy( int level);
    void setLogger( Logger* pLog ) { _pLog = pLog; }
    
private:
    
    Logger *_pLog;
    
    void SQRTOverlap(ub::symmetric_matrix<double> &S, 
                     ub::matrix<double> &Sm2);
    
    bool ProjectExcitons(const ub::matrix<float>& _kap,const ub::matrix<float>& _kbp, 
                         const ub::matrix<float>& ctAB,const ub::matrix<float>& ctBA, 
                         const ub::matrix<float>& _bseA,const ub::matrix<float>& _bseB, 
                         const ub::matrix<float>& _H, ub::matrix<float>& _J );
    
    ub::matrix<float> JAB_singlet;
    ub::matrix<float> JAB_triplet;
    string _spintype;
    bool _doTriplets;
    bool _doSinglets;
    int _levA;
    int _levB;
    int _occA;
    int _unoccA;
    int _occB;
    int _unoccB;
    double      _degeneracy;
    
    
};

}}

#endif	/* _VOTCA_XTP_BSECOUPLING_H */


