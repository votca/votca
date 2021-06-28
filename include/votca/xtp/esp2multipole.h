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
#ifndef VOTCA_XTP_ESP2MULTIPOLE_H
#define VOTCA_XTP_ESP2MULTIPOLE_H

// Standard includes
#include <cstdio>

// Third party includes
#include <boost/filesystem.hpp>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "classicalsegment.h"
#include "espfit.h"
#include "logger.h"
#include "orbitals.h"

namespace votca {
namespace xtp {

class Esp2multipole {
 public:
  Esp2multipole(Logger& log) : log_(log) {
    pairconstraint_.resize(0);
    regionconstraint_.resize(0);
  }

  std::string Identify() { return "esp2multipole"; }

  void Initialize(tools::Property& options);

  StaticSegment Extractingcharges(const Orbitals& orbitals) const;

  std::string GetStateString() const { return state_.ToString(); }

 private:
  void PrintDipoles(const Orbitals& orbitals, const StaticSegment& seg) const;

  QMState state_;
  std::string method_;
  std::string gridsize_;
  bool use_mulliken_;
  bool use_lowdin_;
  bool use_CHELPG_;
  bool do_svd_;
  double conditionnumber_;

  Logger& log_;
  std::vector<std::pair<Index, Index> > pairconstraint_;  //  pairconstraint[i]
                                                          //  is all the
                                                          //  atomindices which
                                                          //  have the same
                                                          //  charge
  std::vector<QMFragment<double> > regionconstraint_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ESP2MULTIPOLE_H
