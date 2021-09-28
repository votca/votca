/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#pragma once
#ifndef VOTCA_XTP_NEIGHBORLIST_H
#define VOTCA_XTP_NEIGHBORLIST_H

// VOTCA includes
#include <votca/tools/globals.h>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/qmnblist.h"
#include "votca/xtp/topology.h"

namespace votca {
namespace xtp {

class Neighborlist final : public QMCalculator {
 public:
  std::string Identify() const { return "neighborlist"; }
  bool WriteToStateFile() const { return true; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Evaluate(Topology& top);

 private:
  Index DetClassicalPairs(Topology& top);

  std::vector<std::string> included_segments_;
  std::map<std::string, std::map<std::string, double> > cutoffs_;
  bool useConstantCutoff_;
  double constantCutoff_;
  bool useExcitonCutoff_;
  double excitonqmCutoff_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_NEIGHBORLIST_H
