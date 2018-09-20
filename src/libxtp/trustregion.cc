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

#include <votca/xtp/trustregion.h>
#include <iostream>

namespace votca {
    namespace xtp {

       Eigen::VectorXd TrustRegion::CalculateStep(const Eigen::VectorXd& gradient, const Eigen::MatrixXd& Hessian, double delta)const{
         //calculate unrestricted step
         Eigen::VectorXd freestep=Hessian.colPivHouseholderQr().solve(-gradient);
         
         const double delta_squared=delta*delta;
         //if inside use the step;
         if(freestep.squaredNorm()<delta_squared){
           return freestep;
         }
         
         //calculate step on the boundary
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hessian);
        const Eigen::VectorXd factor =(es.eigenvectors().transpose()*gradient).cwiseAbs2();
        
        double lambda_s=-es.eigenvalues()(0);
        const double g=gradient.norm();
        const double B1=Hessian.cwiseAbs().maxCoeff();
        double lambda_l=std::max(0.0,std::max(lambda_s,g/delta_squared-B1));
        double lambda_u=g/delta_squared+B1;
         
          // start value for lambda  a bit lower than lowest eigenvalue of Hessian
        double lambda= -es.eigenvalues()(0)+std::sqrt(factor(0))/delta;
        TrustRegionFunction TRF=TrustRegionFunction(factor,es,delta_squared);
        double func_value=0;
        do {
          
           lambda=std::max(lambda,lambda_l);
             lambda=std::min(lambda,lambda_u);
             if(lambda<=lambda_s){
               lambda=std::max(0.001*lambda_u,std::sqrt(lambda_l*lambda_u));
             }
             std::pair<double,double> result=TRF.Evaluate(lambda);
             func_value=result.first;
             
             double update=result.second;
             lambda+=update;
             if(lambda>lambda_s && func_value<0){
               lambda_u=std::min(lambda_u,lambda);
             }else{
               lambda_l=std::max(lambda_l,lambda);
             }
             
             lambda_l=std::max(lambda_l,lambda_s);
            
            }while (std::abs(func_value)>(delta_squared*1e-5));
                
            
            //this is effectively the solution of (H+I\lambda)*\Delta p=-g with \lambda adjusted so that ||p||=delta 
            Eigen::VectorXd new_delta_pos = Eigen::VectorXd::Zero(gradient.size());
            for (int i = 0; i < gradient.size(); i++) {
                new_delta_pos -= es.eigenvectors().col(i) * (es.eigenvectors().col(i).transpose()*gradient) / (es.eigenvalues()(i) + lambda);
            }
            return new_delta_pos;
        }
         
       }
      

     
}
