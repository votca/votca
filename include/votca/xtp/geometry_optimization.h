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

#ifndef __XTP_GEOMETRY_OPTIMIZATION__H
#define __XTP_GEOMETRY_OPTIMIZATION__H


#include <votca/xtp/qmatom.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/segment.h>
#include <stdio.h>
#include <votca/xtp/gwbseengine.h>



namespace votca {
    namespace xtp {



        class GeometryOptimization {
        public:

            GeometryOptimization(GWBSEENGINE& gwbse_engine, QMPackage *qmpackage, 
                    std::vector<ctp::Segment*> segments, Orbitals& orbitals) : 
            _gwbse_engine(gwbse_engine), _qmpackage(qmpackage), _segments(segments), _orbitals(orbitals) {
                
            };

            ~GeometryOptimization() {
            };

            void BFGSStep(int& _iteration, bool& _update_hessian, Eigen::MatrixXd& _force, 
            Eigen::MatrixXd& _force_old, Eigen::MatrixXd& _current_xyz, 
            Eigen::MatrixXd& _old_xyz, Eigen::MatrixXd& _hessian, Eigen::MatrixXd& _xyz_shift, Eigen::MatrixXd& _trial_xyz);
            void Initialize(tools::Property *options);

            void setLog(ctp::Logger* pLog) {
                _pLog = pLog;
            }

            void Evaluate();



        private:

            int _natoms;
            int _nsegments;


            int _opt_state;
            std::string _spintype;
            std::string _forces;
            std::string _opt_type;
            std::string _optimizer;
            std::string _force_method;

            GWBSEENGINE _gwbse_engine;
            QMPackage* _qmpackage;
            std::vector<ctp::Segment*> _segments;
            Orbitals& _orbitals;

            tools::Property _optimizer_options;
            tools::Property _force_options;

            ctp::Logger *_pLog;
        };

    }
}
#endif /* GEOMETRY_OPTIMIZATION_H */
