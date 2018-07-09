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

#ifndef _VOTCA_XTP_PPM_H
#define _VOTCA_XTP_PPM_H
#include <votca/xtp/eigen.h>
#include <votca/xtp/rpa.h>


namespace votca {
namespace xtp {




class PPM {
 public:
     
     PPM():screening_r(0.0),screening_i(0.5){};
     
    //This sets the screening frequencies for real and imaginary part in hartree
     

double getScreening_r()const{return screening_r;}

double getScreening_i()const{return screening_i;}
 void PPM_construct_parameters(const RPA& rpa);

 const Eigen::VectorXd& getPpm_weight() const {
     return _ppm_weight;
 }

 const Eigen::VectorXd& getPpm_freq() const {
     return _ppm_freq;
 }

 const Eigen::MatrixXd& getPpm_phi() const {
     return _ppm_phi;
 }     
     
 void FreeMatrix(){
     _ppm_phi.resize(0,0);
 }
     
 private:
 
     double screening_r;
     double screening_i;

  // PPM related variables and functions
  Eigen::MatrixXd _ppm_phi;
  Eigen::VectorXd _ppm_freq;
  Eigen::VectorXd _ppm_weight;
  
  
 

  

 
};
}
}

#endif /* _VOTCA_XTP_GWBSE_H */
