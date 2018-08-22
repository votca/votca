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

            BFGSTRM(GWBSEEngine& gwbse_engine, QMPackage* qmpackage, Orbitals& orbitals, Forces& force_engine)
            : _gwbse_engine(gwbse_engine), _qmpackage(qmpackage), _orbitals(orbitals), _force_engine(force_engine), _iteration(0) {
            };

            ~BFGSTRM() {
            };

            int Iteration() {
                return _iteration;
            };
            void Initialize(tools::Property &options);
            void Checkpoint(std::vector<ctp::Segment* >& _molecule);
            void WriteIteration(FILE* out, ctp::Segment* _segment);

            void setLog(ctp::Logger* pLog) {
                _pLog = pLog;
            }

            void Optimize();



        private:
            
            GWBSEEngine& _gwbse_engine;
            QMPackage* _qmpackage;
            Orbitals& _orbitals;
            Forces _force_engine;

            unsigned _natoms;
            unsigned _iteration;
            Eigen::MatrixX3d _force;
            Eigen::MatrixX3d _force_old;
            Eigen::MatrixX3d _xyz_shift;
            Eigen::MatrixX3d _current_xyz;
            Eigen::MatrixX3d _old_xyz;
            Eigen::MatrixX3d _trial_xyz;
            Eigen::MatrixXd _hessian;

            bool _step_accepted;
            bool _update_hessian;
            bool _restart_opt;

            QMState _opt_state;
            double _displacement;
            double _convergence;
            double _RMSForce_convergence;
            double _MaxForce_convergence;
            double _RMSStep_convergence;
            double _MaxStep_convergence;
            double _trust_radius;
            double _trust_radius_max;
            double _delta_energy_estimate;
            double _norm_delta_pos;
            std::string _forces;
            std::string _opt_type;
            std::string _optimizer;
            std::string _force_method;
            unsigned _max_iteration;

            

            tools::Property _optimizer_options;
            tools::Property _force_options;

            ctp::Logger *_pLog;

            void BFGSStep();
            void Rewrite2Vectors();
            void Rewrite2Matrices();
            void UpdateHessian();
            void PredictDisplacement();
            void RegularizeStep();
            void QuadraticEnergy();
            void UpdateCoordinatesOrbitals();
            void Report();
            void OrbitalsToMatrix();
            void AcceptReject();
            void WriteTrajectory();

            double GetEnergy();

            bool OutsideTrustRegion(const double& _step);
            bool GeometryConverged();

            std::string Converged(bool converged);


            // vector storage for steps, let's rethink that later
            unsigned _dim;
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
