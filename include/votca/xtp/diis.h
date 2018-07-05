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

#ifndef _VOTCA_XTP_DIIS__H
#define _VOTCA_XTP_DIIS__H



#include <votca/xtp/eigen.h>
#include <vector>

namespace votca { namespace xtp {
 
 
 class DIIS{
public:

    DIIS( ):success(true) { };
   ~DIIS() {
       
    for (std::vector< Eigen::MatrixXd* >::iterator it = _errormatrixhist.begin() ; it !=_errormatrixhist.end(); ++it){
         delete *it;
     }
    _errormatrixhist.clear();    
       
   for (std::vector< std::vector<double>* >::iterator it = _Diis_Bs.begin() ; it !=_Diis_Bs.end(); ++it){
         delete *it;
     }
    _Diis_Bs.clear(); 
   };
   
    void Update(unsigned _maxerrorindex, const Eigen::MatrixXd& errormatrix);
    Eigen::VectorXd CalcCoeff();
    
    void setHistLength(unsigned length){_histlength=length;}
   
    bool Info(){return success;}
    
 private:
     
     
     
     
     bool success;
     unsigned _histlength;
     std::vector< std::vector<double>* >  _Diis_Bs;
     std::vector< Eigen::MatrixXd* >   _errormatrixhist;
  
 };
 
 
    
}}

#endif	

