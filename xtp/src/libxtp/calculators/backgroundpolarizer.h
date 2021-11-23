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
#ifndef VOTCA_XTP_BACKGROUNDPOLARIZER_H
#define VOTCA_XTP_BACKGROUNDPOLARIZER_H

// Local VOTCA includes
#include "votca/xtp/background.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/topology.h"

namespace votca {
namespace xtp {


/*!
 * \brief class that builds and polarizes the periodic background using an Ewald
 * calculation.
 */
class BackgroundPolarizer final : public QMCalculator {
 public:
  BackgroundPolarizer() = default;
  ~BackgroundPolarizer() = default;

  std::string Identify() const { return "backgroundpolarizer"; }

  bool WriteToStateFile() const { return true; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Evaluate(Topology &top);

 private:
  Logger log_;
  EwaldOptions ewd_options_;
  std::string mapfile_;
  std::string output_file_name_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDBG_H
