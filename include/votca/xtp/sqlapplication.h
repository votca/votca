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

#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/topology.h>

#include <votca/xtp/statesaversqlite.h>
#include <votca/xtp/xtpapplication.h>

namespace votca {
namespace xtp {

class SqlApplication : public XtpApplication {
 public:
  SqlApplication();

  ~SqlApplication(){};

  void Initialize();
  bool EvaluateOptions();
  void Run(void);

  virtual void BeginEvaluate(int nThreads);
  virtual bool EvaluateFrame();
  virtual void EndEvaluate();

  void AddCalculator(QMCalculator *calculator);

 protected:
  std::vector<std::unique_ptr<QMCalculator> > _calculators;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SQLAPPLICATION_H
