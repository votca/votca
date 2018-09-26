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

#ifndef XTP_GEOMETRY_OPTIMIZATION_H
#define XTP_GEOMETRY_OPTIMIZATION_H


#include <votca/xtp/qmatom.h>
#include <votca/xtp/logger.h>
#include <stdio.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/bfgs-trm.h>
#include <votca/xtp/energy_costfunction.h>


namespace votca {
    namespace xtp {

        class GeometryOptimization {
        public:

            GeometryOptimization(GWBSEEngine& gwbse_engine,Orbitals& orbitals) :
            _gwbse_engine(gwbse_engine), _orbitals(orbitals) {

            };

            void Initialize(tools::Property& options);

            void setLog(Logger* pLog) {
                _pLog = pLog;
            }

            void Evaluate();

        private:

            static void Report(const BFGSTRM& bfgstrm, Logger* pLog);
            static void WriteTrajectory(const std::string& filename, std::vector< QMAtom* >& atoms,
                                        const BFGSTRM& bfgstrm);

            QMState _opt_state;
            std::string _optimizer;
            std::string _trajfile;
            GWBSEEngine& _gwbse_engine;
            Orbitals& _orbitals;

            Energy_costfunction::conv_paras _conv;
            int _max_iteration;
            double _trust_radius;

            tools::Property _filter_options;
            tools::Property _force_options;

            Logger *_pLog;
        };

    }
}
#endif // VOTCA_XTP_GEOMETRY_OPTIMIZATION_H 
