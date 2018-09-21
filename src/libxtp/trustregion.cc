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
#include <eigen3/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h>
#include <iomanip>

namespace votca {
    namespace xtp {

       Eigen::VectorXd TrustRegion::CalculateStep(const Eigen::VectorXd& gradient, const Eigen::MatrixXd& Hessian, double delta)const{
         //calculate unrestricted step
         Eigen::VectorXd freestep=Hessian.colPivHouseholderQr().solve(-gradient);

         //if inside use the step;
         if(freestep.norm()<delta){
           return freestep;
         }
         
         //calculate step on the boundary
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hessian);
        const Eigen::VectorXd factor =(es.eigenvectors().transpose()*gradient).cwiseAbs2();
        double lambda=0;
        //hard case
        if(std::abs(factor[0])<0){
            lambda=-es.eigenvalues()(0);
            int size=factor.size()-1;
            Eigen::ArrayXd factor_small=factor.tail(size).array();
            Eigen::ArrayXd quotient=es.eigenvalues().tail(size).array()+lambda;
            const double p2=(factor_small/(quotient.pow(2))).sum();
            const double tau=std::sqrt(delta*delta-p2);
            Eigen::VectorXd new_delta_pos = -tau*es.eigenvectors().col(0);
             for (int i = 1; i < gradient.size(); i++) {
                new_delta_pos -= es.eigenvectors().col(i) * (es.eigenvectors().col(i).transpose()*gradient) / (es.eigenvalues()(i) + lambda);
            }
            
            return new_delta_pos; 
        }

        
          // start value for lambda  a bit higher than lowest eigenvalue of Hessian
        lambda= -es.eigenvalues()(0)+std::sqrt(factor(0))/delta;
        TrustRegionFunction TRF=TrustRegionFunction(factor,es,delta);
        double func_value=0;
        do {
            std::pair<double,double> result=TRF.Evaluate(lambda);
             func_value=result.first;
             double update=result.second;
             lambda+=update;
            }while (std::abs(func_value)>(1/delta*1e-9));
            
            //this is effectively the solution of (H+I\lambda)*\Delta p=-g with \lambda adjusted so that ||p||=delta 
            Eigen::VectorXd new_delta_pos = Eigen::VectorXd::Zero(gradient.size());
            for (int i = 0; i < gradient.size(); i++) {
                new_delta_pos -= es.eigenvectors().col(i) * (es.eigenvectors().col(i).transpose()*gradient) / (es.eigenvalues()(i) + lambda);
            }
            return new_delta_pos;
        }
         
    } 
}
