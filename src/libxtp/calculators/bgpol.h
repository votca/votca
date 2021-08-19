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
#ifndef VOTCA_XTP_BGPOL_H
#define VOTCA_XTP_BGPOL_H

// Local VOTCA includes
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/classicalsegment.h"
#include "ewald/backgroundpolarizer.h"

namespace votca {
namespace xtp {

class BGPol final : public QMCalculator {
 public:
  BGPol() = default;
  ~BGPol() = default;

  std::string Identify() const { return "bgpol"; }

  bool WriteToStateFile() const { return true; }

  void WriteToHdf5(std::string filename, std::vector<EwdSegment> background) const;

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Evaluate(Topology &top);

 private:
 EwaldOptions ewd_options;
 std::string _mapfile;
 Logger _log;
 std::vector<PolarSegment> _polar_background;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDBG_H
