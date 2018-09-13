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

#ifndef VOTCA_XTP_ADIIS_COSTFUNCTION_H
#define VOTCA_XTP_ADIIS_COSTFUNCTION_H


#include <votca/xtp/basisset.h>
#include "ceres/ceres.h"


namespace votca {
    namespace xtp {

        class ADIIS_costfunction : public ceres::FirstOrderFunction {
        public:

            ADIIS_costfunction(Eigen::VectorXd DiF, Eigen::MatrixXd DiFj) {
                _DiF = DiF;
                _DiFj = DiFj;
            }

            virtual bool Evaluate(const double* parameters,
                    double* cost,
                    double* gradient) const {
                Eigen::Map<const Eigen::VectorXd> x(parameters, _DiF.size());
                Eigen::VectorXd c = x.cwiseAbs2();
                double xnorm = c.sum();
                c /= xnorm;

                cost[0] = (2 * c.transpose() * _DiF + c.transpose() * _DiFj * c).value();
                if (gradient != NULL) {
                    Eigen::VectorXd dEdc = 2.0 * _DiF + _DiFj * c + _DiFj.transpose() * c;
                    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(c.size(), c.size());
                    for (int i = 0; i < jac.rows(); i++) {
                        for (int j = 0; j < jac.cols(); j++) {
                            jac(i, j) = -c(i)*2.0 * x(j) / xnorm;
                        }
                        // Extra term on diagonal
                        jac(i, i) += 2.0 * x(i) / xnorm;
                    }
                    Eigen::VectorXd dEdxv = jac.transpose() * dEdc;
                    for (int i = 0; i < dEdxv.size(); ++i) {
                        gradient[i] = dEdxv(i);
                    }
                }
                return true;
            }

            virtual int NumParameters() const {
                return _DiF.size();
            }


        private:
            Eigen::VectorXd _DiF;
            Eigen::MatrixXd _DiFj;


        };

    }
}
#endif // VOTCA_XTP_ADIIS_COSTFUNCTION_H
