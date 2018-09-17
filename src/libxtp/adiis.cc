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
#include "votca/xtp/adiis.h"
#include "ceres/ceres.h"
#include <votca/xtp/adiis_costfunction.h>


namespace votca { namespace xtp {
  
      
      
   Eigen::VectorXd ADIIS::CalcCoeff(const std::vector< Eigen::MatrixXd* >& dmathist,const std::vector< Eigen::MatrixXd* >& mathist){
      success=true;
      int size=dmathist.size();
      
      const Eigen::MatrixXd& dmat=*dmathist[size-1];
      const Eigen::MatrixXd& H=*mathist[size-1];
      Eigen::VectorXd DiF = Eigen::VectorXd::Zero(size);
      Eigen::MatrixXd DiFj = Eigen::MatrixXd::Zero(size, size);


      for (int i = 0; i < size; i++) {
        DiF(i) = ((*dmathist[i]) - dmat).cwiseProduct(H).sum();
      }

      for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
          DiFj(i, j) = ((*dmathist[i]) - dmat).cwiseProduct((*mathist[j]) - H).sum();
        }
      }   
     
 
   ceres::GradientProblem problem(new ADIIS_costfunction(DiF,DiFj));
   // Starting point: equal weights on all matrices
   Eigen::VectorXd coeffs=Eigen::VectorXd::Constant(size,1.0/size);

   ceres::GradientProblemSolver::Options options;
   options.minimizer_progress_to_stdout=false;
   options.logging_type=ceres::LoggingType::SILENT;
   options.gradient_tolerance=1e-8;
   options.max_num_iterations=1000;
   ceres::GradientProblemSolver::Summary summary;
   ceres::Solve(options,problem,coeffs.data(),&summary);
   //std::cout << summary.FullReport() << "\n";
  success=summary.IsSolutionUsable();
    
  coeffs=coeffs.cwiseAbs2();
  double xnorm=coeffs.sum();
  coeffs/=xnorm;

  
  if(std::abs(coeffs.tail(1).value())<0.001){     
        success=false;
      }
  return coeffs;
}

}}
