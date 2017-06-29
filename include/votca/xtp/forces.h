/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/ctp/qmatom.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/segment.h>
#include <stdio.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/qminterface.h>



using namespace std;

namespace votca {
    namespace xtp {

        namespace ub = boost::numeric::ublas;

        class Forces {
        public:

            Forces(GWBSEENGINE& gwbse_engine, QMPackage* qmpackage, vector<ctp::Segment*> segments, Orbitals* orbitals)
            : _gwbse_engine(gwbse_engine), _qmpackage(qmpackage), _segments(segments), _orbitals(orbitals), _remove_total_force(false), _remove_CoM_force(false) {
            };

            ~Forces() {
            };

            void Initialize(Property *options);
            void Calculate(const double& energy);

            void NumForceForward(double energy, std::vector< ctp::Atom* > ::iterator ait, ub::matrix_range< ub::matrix<double> >& _force,
                    std::vector<ctp::Segment*> _molecule);
            void NumForceCentral(double energy, std::vector< ctp::Atom* > ::iterator ait, ub::matrix_range< ub::matrix<double> >& _force,
                    std::vector<ctp::Segment*> _molecule);

            void setLog(ctp::Logger* pLog) {
                _pLog = pLog;
            }

            void SetSpinType(const string spin_type) {
                _spin_type = spin_type;
            };

            void SetOptState(const int opt_state) {
                _opt_state = opt_state;
            };

            string GetSpinType() {
                return _spin_type;
            };

            int GetOptState() {
                return _opt_state;
            };

            ub::matrix<double> GetForces() {
                return _forces;
            };
            void Report();




        private:

            double _displacement;
            string _force_method;
            string _spin_type;

            bool _remove_total_force;
            bool _remove_CoM_force;

            int _nsegments;
            int _natoms;
            int _opt_state;

            GWBSEENGINE _gwbse_engine;
            QMPackage* _qmpackage;
            vector<ctp::Segment*> _segments;
            Orbitals* _orbitals;

            ub::matrix<double> _forces;

            Property _force_options;

            void RemoveTotalForce();
            void RemoveCoMForce();
            ub::vector<double> TotalForce();

            QMMInterface _qminterface;
            ctp::Logger *_pLog;
        };

    }
}
#endif /* FORCES_H */
