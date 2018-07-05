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

#ifndef _VOTCA_XTP_MIXING__H
#define _VOTCA_XTP_MIXING__H


#include <votca/xtp/aomatrix.h>
#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>



namespace votca { namespace xtp {

  
 //Mixing according to Zerner, M.C., Hehenberger, M., 1979. 
 //A dynamical damping scheme for converging molecular scf calculations. 
 //Chemical Physics Letters 62, 550–554. https://doi.org/10.1016/0009-2614(79)80761-7

 class Mixing{
public:

    Mixing() {};
    
    
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
   
   void Configure(double mixingparameter,const  Eigen::MatrixXd* S){
       _S=S;
        _mixingparameter=mixingparameter;
        if(_mixingparameter>0 && _mixingparameter<1.0){
            _automaticmixing=false;
        }else{
            _automaticmixing=true;
        }
   }
   
   double getAlpha(){return _alpha;}
   
  Eigen::MatrixXd MixDmat(const Eigen::MatrixXd& dmatin,const Eigen::MatrixXd& dmatout);
   void Updatemix(const Eigen::MatrixXd& dmatin,const Eigen::MatrixXd& dmatout ); 
   
 private:
     
     
    Eigen::VectorXd Mullikencharges(const Eigen::MatrixXd& dmat);
    const Eigen::MatrixXd*_S;
    double _alpha;
    bool _automaticmixing;
    double _mixingparameter;
    std::vector< Eigen::VectorXd* >      _Pin;
    std::vector< Eigen::VectorXd* >      _Pout;
    
  
 };
    
}}

#endif	

