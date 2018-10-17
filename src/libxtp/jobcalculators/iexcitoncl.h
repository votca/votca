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

#ifndef _CALC_COUPLING_EXCL_H
#define _CALC_COUPLING_EXCL_H

#include <votca/tools/property.h>

#include <votca/xtp/parallelxjobcalc.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <votca/xtp/qmstate.h>

namespace votca {
    namespace xtp {

        /**
         * \brief Evaluates Transition Charge distributions classically
         *
         * Evaluates the electrostatic classical coupling between molecules in 
         * their excited states.
         * Callname: iexcitoncl
         */

        class IEXCITON : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult > {
        public:

            void Initialize(tools::Property *options);

            string Identify() {
                return "iexcitoncl";
            }

            Job::JobResult EvalJob(Topology *top, Job *job, QMThread *Thread);

            void WriteJobFile(Topology *top);
            void ReadJobFile(Topology *top);

        private:
            QMState GetElementFromMap(const std::string& elementname )const;
            std::map<std::string, QMState> FillParseMaps(const string& Mapstring);
            double _cutoff;
            bool _induce;
            std::map<std::string,QMState> _statemap;
            string _emp_file;
            string _xml_file;
            void PreProcess(Topology *top);
            

        };

    }
}
#endif /* _CALC_INTEGRALS_DFT_H */
