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

#ifndef __XTP_ADIIS_COSTFUNCTION__H
#define __XTP_ADIIS_COSTFUNCTION__H


#include <votca/xtp/basisset.h>


namespace votca {
    namespace xtp {

        class ADIIS_costfunction :  public Optimiser_costfunction {
        public:

            ADIIS_costfunction(Eigen::VectorXd DiF, Eigen::MatrixXd DiFj) {
                _DiF = DiF;
                _DiFj = DiFj;
            }
            
             double EvaluateCost(const Eigen::VectorXd& parameters){
                 Eigen::VectorXd c = parameters.cwiseAbs2();
                double xnorm = c.sum();
                c /= xnorm;
                return (2 * c.transpose() * _DiF + c.transpose() * _DiFj * c).value();
             }
            
            Eigen::VectorXd EvaluateGradient(const Eigen::VectorXd& parameters){
                Eigen::VectorXd c = parameters.cwiseAbs2();
                double xnorm = c.sum();
                c /= xnorm;
                Eigen::VectorXd dEdc = 2.0 * _DiF + _DiFj * c + _DiFj.transpose() * c;
                Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(c.size(), c.size());
                for (int i = 0; i < jac.rows(); i++) {
                    for (int j = 0; j < jac.cols(); j++) {
                        jac(i, j) = -c(i)*2.0 * parameters(j) / xnorm;
                    }
                    // Extra term on diagonal
                    jac(i, i) += 2.0 * parameters(i) / xnorm;
                }
                return jac.transpose() * dEdc; 
            }

            int NumParameters() const {
                return _DiF.size();
            }
            
            bool Converged(const Eigen::VectorXd& delta_parameters,
                    double delta_cost, const Eigen::VectorXd& gradient){
                return gradient.cwiseAbs().sum()<1.e-9;
            }


        private:
            Eigen::VectorXd _DiF;
            Eigen::MatrixXd _DiFj;


        };

    }
}
#endif /* FORCES_H */
