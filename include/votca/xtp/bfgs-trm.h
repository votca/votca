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

#ifndef __XTP_BFGSTRM__H
#define __XTP_BFGSTRM__H


#include <votca/ctp/logger.h>
#include <votca/ctp/segment.h>
#include <stdio.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/forces.h>






namespace votca {
    namespace xtp {

       
        class BFGSTRM {
        public:

            BFGSTRM(GWBSEEngine& gwbse_engine, Statefilter& filter, Orbitals& orbitals, Forces& force_engine)
            : _gwbse_engine(gwbse_engine), _filter(filter), _orbitals(orbitals), _force_engine(force_engine), _iteration(0) {
            };

            int Iteration() const{
                return _iteration;
            };
            void Initialize(tools::Property &options);

            void setLog(ctp::Logger* pLog) {
                _pLog = pLog;
            }

            void Optimize();

        private:
            
            void BFGSStep();
            void WriteMatrixToVector();
            void WriteCoordinates2Matrices();
            void UpdateHessian();
            Eigen::VectorXd PredictDisplacement();
            void RegularizeStep();
            double QuadraticEnergy();
            void Vector2QMAtoms();
            void Report();
            Eigen::VectorXd QMAtomsToVector(std::vector<QMAtoms*>& vec);
            Eigen::VectorXd Write3XMatrixToVector(const Eigen::MatrixX3d& matrix);
            bool AcceptRejectStep();
            void WriteTrajectory();

            double GetEnergy();

            bool OutsideTrustRegion(double step);
            bool CheckConvergence();

            std::string Converged(bool converged);
            
            GWBSEEngine& _gwbse_engine;
            Statefilter& _filter;
            Orbitals& _orbitals;
            Forces& _force_engine;

            unsigned _iteration;
            Eigen::MatrixXd _hessian;

            bool _update_hessian;

            double _displacement;
            double _convergence;
            double _RMSForce_convergence;
            double _MaxForce_convergence;
            double _RMSStep_convergence;
            double _MaxStep_convergence;
            double _trust_radius;

            unsigned _max_iteration;

            ctp::Logger *_pLog;

            Eigen::VectorXd _previous_pos;
            Eigen::VectorXd _current_pos;
            Eigen::VectorXd _previous_gradient;
            Eigen::VectorXd _current_gradient;
            Eigen::VectorXd _delta_pos;
            double _new_energy;
            double _last_energy;
            double _energy_delta;

            // convergence
            bool _energy_converged;
            bool _RMSForce_converged;
            bool _MaxForce_converged;
            bool _RMSStep_converged;
            bool _MaxStep_converged;
            double _RMSForce;
            double _MaxForce;
            double _RMSStep;
            double _MaxStep;


        };

    }
}
#endif /* BFGSTRM_H */
