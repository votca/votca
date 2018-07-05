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

#ifndef _VOTCA_XTP_ADIIS__H
#define _VOTCA_XTP_ADIIS__H



#include <votca/xtp/eigen.h> 
#include <vector>

namespace votca { namespace xtp {


 
 class ADIIS{
public:

    ADIIS():success(true) {};
   ~ADIIS() {};
   
    Eigen::VectorXd CalcCoeff(const std::vector< Eigen::MatrixXd* >& _dmathist,const std::vector< Eigen::MatrixXd* >& _mathist);
    
   
   bool Info(){return success;}
 private:
     
     bool success;
  
 };
    
}}

#endif	

