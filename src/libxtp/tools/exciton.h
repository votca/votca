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

#ifndef _VOTCA_XTP_EXCITON_H
#define _VOTCA_XTP_EXCITON_H

#include <stdio.h>

#include <votca/ctp/logger.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/qmtool.h>
#include <votca/ctp/segment.h>

#include <votca/tools/linalg.h>
#include <votca/tools/constants.h>


namespace votca {
    namespace xtp {
        using namespace std;

        class Exciton : public ctp::QMTool {
        public:

            Exciton() {
            };

            ~Exciton() {
            };

            string Identify() {
                return "exciton";
            }

            void Initialize(Property *options);
            bool Evaluate();




        private:


            string _xyzfile;
            string _xml_output;    // .xml output
            string _package;
            string _archive_file; // .orb file to parse to
            string _reporting;
            string _guess_orbA;
            string _guess_orbB;

            Property _package_options;
            Property _gwbseengine_options;
            Property _geoopt_options;

            ctp::Logger _log;

            bool _do_optimize;

            void ReadXYZ(ctp::Segment* _segment, string filename);

        };



    }
}


#endif
