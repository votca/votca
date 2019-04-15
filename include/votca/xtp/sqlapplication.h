/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef VOTCA_XTP_SQLAPPLICATION_H
#define VOTCA_XTP_SQLAPPLICATION_H

#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/topology.h>

#include <votca/xtp/statesaversqlite.h>
#include <votca/xtp/xtpapplication.h>

namespace votca {
namespace xtp {

class SqlApplication : public XtpApplication {
 public:
  SqlApplication();

  ~SqlApplication() {
    for (ctp::QMCalculator* calculator : _calculators) {
      delete calculator;
    }
  };

  void Initialize();
  bool EvaluateOptions();
  void Run(void);

  virtual void BeginEvaluate(int nThreads);
  virtual bool EvaluateFrame();
  virtual void EndEvaluate();

  void AddCalculator(ctp::QMCalculator* calculator);

 protected:
  ctp::Topology _top;
  list<ctp::QMCalculator*> _calculators;
};

}  // namespace xtp
}  // namespace votca

#endif /* _QMApplication_H */
