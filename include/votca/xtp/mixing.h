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

#ifndef _VOTCA_XTP_MIXING__H
#define _VOTCA_XTP_MIXING__H

#include <votca/tools/linalg.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>



namespace votca { namespace xtp {
 namespace ub = boost::numeric::ublas;
  
 
 class Mixing{
public:

    Mixing(bool automaticmixing,double mixingparameter,Eigen::MatrixXd* _S,votca::ctp::Logger *pLog) {
        _mixingparameter=mixingparameter;
        _automaticmixing=automaticmixing;
        S=_S;
        _pLog=pLog;
    };
   ~Mixing() {
    for (std::vector< Eigen::VectorXd* >::iterator it = _Pout.begin() ; it !=_Pout.end(); ++it){
         delete *it;
     }
     _Pout.clear();
     for (std::vector< Eigen::VectorXd* >::iterator it = _Pin.begin() ; it !=_Pin.end(); ++it){
         delete *it;
     }
     _Pin.clear();
   }
   
  
   
  Eigen::MatrixXd MixDmat(const Eigen::MatrixXd& dmatin,const Eigen::MatrixXd& dmatout,bool noisy=true );
   void Updatemix(const Eigen::MatrixXd& dmatin,const Eigen::MatrixXd& dmatout ); 
   
 private:
     
     
    
   Eigen::VectorXd Mullikencharges(const Eigen::MatrixXd& dmat);
  
    ctp::Logger *_pLog;
    Eigen::MatrixXd* S;
    bool _automaticmixing;
    double _mixingparameter;
    std::vector< Eigen::VectorXd* >      _Pin;
    std::vector< Eigen::VectorXd* >      _Pout;
    
  
 };
    
}}

#endif	

