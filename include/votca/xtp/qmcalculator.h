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
/// For an earlier history see ctp repo commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_QMCALCULATOR_H
#define VOTCA_XTP_QMCALCULATOR_H

#include <votca/tools/calculator.h>

namespace TOOLS = votca::tools;

namespace votca { namespace xtp {

class Topology;

class QMCalculator : public TOOLS::Calculator{
public:

                    QMCalculator() {}
    virtual        ~QMCalculator() {}

    virtual std::string  Identify() = 0;

    virtual void    Initialize(TOOLS::Property *options) = 0;
    virtual bool    EvaluateFrame(Topology *top) { return true; }
    virtual void    EndEvaluate(Topology *top) { }

};

}}

#endif // VOTCA_XTP_QMCALCULATOR_H 
