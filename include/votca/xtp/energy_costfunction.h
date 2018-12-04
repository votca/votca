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

#ifndef __XTP_ENERGY_COSTFUNCTION__H
#define __XTP_ENERGY_COSTFUNCTION__H

#include <votca/xtp/optimiser_costfunction.h>

#include "orbitals.h"
#include "forces.h"

namespace votca {
    namespace xtp {

        class Energy_costfunction : public Optimiser_costfunction {
        public:

            struct conv_paras {
                double deltaE;
                double RMSForce;
                double MaxForce;
                double RMSStep;
                double MaxStep;
                int maxforceindex=0;
                int maxstepindex=0;
            };

            Energy_costfunction(GWBSEEngine& gwbse_engine, Statefilter& filter,
                    Orbitals& orbitals, Forces& force_engine)
            : _gwbse_engine(gwbse_engine), _filter(filter),
            _orbitals(orbitals), _force_engine(force_engine) {
            };

            double EvaluateCost(const Eigen::VectorXd& parameters);
            
            Eigen::VectorXd EvaluateGradient(const Eigen::VectorXd& parameters);

            int NumParameters() const {
                return _orbitals.QMAtoms().size()*3;
            };

            bool Converged(const Eigen::VectorXd& delta_parameters,
                    double delta_cost, const Eigen::VectorXd& gradient);

            void ForcesReport()const {
                return _force_engine.Report();
            }
            
            const conv_paras& getConvParas()const {
                return _convpara;
            }

            void setConvergenceParameters(const conv_paras& convergence) {
                _convpara = convergence;
            }
            
            void setLog(ctp::Logger* pLog) {
                _pLog = pLog;
            }

            
            void Report(const conv_paras& val);
            static void Vector2QMAtoms(const Eigen::VectorXd& pos, std::vector<QMAtom*>& atoms);
            static Eigen::VectorXd QMAtoms2Vector(std::vector<QMAtom*>& atoms);
            static Eigen::VectorXd Write3XMatrixToVector(const Eigen::MatrixX3d& matrix);

        private:
            
            static std::string Converged(double val, double limit);
            GWBSEEngine& _gwbse_engine;
            Statefilter& _filter;
            Orbitals& _orbitals;
            Forces& _force_engine;
            int _iteration=0;
            double _energy;
            
            conv_paras _convpara;

            
            ctp::Logger *_pLog;

            
            


        };

    }
}
#endif /* FORCES_H */
