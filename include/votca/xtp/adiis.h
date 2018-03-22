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

#ifndef _VOTCA_XTP_ADIIS__H
#define _VOTCA_XTP_ADIIS__H



#include <votca/xtp/basisset.h> 
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>







namespace votca { namespace xtp {


 
 class ADIIS{
public:

    ADIIS():success(true) {};
   ~ADIIS() {};
   
    Eigen::VectorXd CalcCoeff(const std::vector< Eigen::MatrixXd* >& _dmathist,const std::vector< Eigen::MatrixXd* >& _mathist);
    
    double get_E_adiis(const gsl_vector * x) const;

    void get_dEdx_adiis(const gsl_vector * x, gsl_vector * dEdx) const;
    void get_E_dEdx_adiis(const gsl_vector * x, double * Eval, gsl_vector * dEdx) const;
   bool Info(){return success;}
 private:
     
     bool success;
    Eigen::VectorXd                  _DiF;
    Eigen::MatrixXd                  _DiFj;
  
   
    
    
 Eigen::VectorXd compute_c(const gsl_vector * x);
 /// Compute jacobian
 Eigen::MatrixXd compute_jac(const gsl_vector * x);
 /// Compute energy
 double min_f(const gsl_vector * x, void * params);
 /// Compute derivative
 void min_df(const gsl_vector * x, void * params, gsl_vector * g);
 /// Compute energy and derivative
void min_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g);
  
 };
 
 
 namespace adiis {
  /// Compute weights
  Eigen::VectorXd compute_c(const gsl_vector * x);
  /// Compute jacobian
  Eigen::MatrixXd compute_jac(const gsl_vector * x);

  /// Compute energy
  double min_f(const gsl_vector * x, void * params);
  /// Compute derivative
  void min_df(const gsl_vector * x, void * params, gsl_vector * g);
  /// Compute energy and derivative
  void min_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g);
};
    
}}

#endif	

