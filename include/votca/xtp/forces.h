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

#ifndef __XTP_FORCES__H
#define __XTP_FORCES__H


#include <votca/xtp/qmatom.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/segment.h>
#include <stdio.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/qminterface.h>



namespace votca {
    namespace xtp {


        class Forces {
        public:

            Forces(GWBSEENGINE& gwbse_engine, QMPackage* qmpackage, std::vector<xtp::Segment*> segments, Orbitals* orbitals)
            : _gwbse_engine(gwbse_engine), _qmpackage(qmpackage), _segments(segments), _orbitals(orbitals), _remove_total_force(false), _remove_CoM_force(false) {
            };

            ~Forces() {
            };

            void Initialize(tools::Property *options);
            void Calculate(const double& energy);

            Eigen::Vector3d NumForceForward(double energy, std::vector< xtp::Atom* > ::iterator ait,std::vector<xtp::Segment*> _molecule);
            Eigen::Vector3d NumForceCentral(double energy, std::vector< xtp::Atom* > ::iterator ait,std::vector<xtp::Segment*> _molecule);

            void setLog(xtp::Logger* pLog) {
                _pLog = pLog;
            }

            void SetSpinType(const std::string spin_type) {
                _spin_type = spin_type;
            };

            void SetOptState(const int opt_state) {
                _opt_state = opt_state;
            };

            std::string GetSpinType() {
                return _spin_type;
            };

            int GetOptState() {
                return _opt_state;
            };

            Eigen::MatrixXd GetForces() {
                return _forces;
            };
            void Report();




        private:

            double _displacement;
            std::string _force_method;
            std::string _spin_type;

            
            
            bool _noisy_output;

            int _nsegments;
            unsigned _natoms;
            int _opt_state;

            GWBSEENGINE _gwbse_engine;
            QMPackage* _qmpackage;
            std::vector<xtp::Segment*> _segments;
            Orbitals* _orbitals;
            bool _remove_total_force;
            bool _remove_CoM_force;

            Eigen::MatrixX3d _forces;

            tools::Property _force_options;

            void RemoveTotalForce();
            void RemoveCoMForce();
            Eigen::Vector3d TotalForce();

            QMMInterface _qminterface;
            xtp::Logger *_pLog;
        };

    }
}
#endif /* FORCES_H */
