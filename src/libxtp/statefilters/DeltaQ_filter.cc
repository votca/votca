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

// Local VOTCA includes
#include "votca/xtp/populationanalysis.h"

// Local private VOTCA includes
#include "DeltaQ_filter.h"

namespace votca {
namespace xtp {

void DeltaQ_filter::Initialize(const tools::Property& options) {
  std::string indices = options.get("fragment").as<std::string>();
  fragment_ = QMFragment<double>(0, indices);
  fragment_.value() = options.get("threshold").as<double>();
}

void DeltaQ_filter::Info(Logger& log) const {
  XTP_LOG(Log::error, log) << "Using Delta Q tracker for fragment " << fragment_
                           << std::flush;
}

void DeltaQ_filter::UpdateHist(const Orbitals&, QMState) { return; }

std::vector<Index> DeltaQ_filter::CalcIndeces(const Orbitals& orb,
                                              QMStateType type) const {

  if (!type.isExciton()) {
    throw std::runtime_error("ChargeTransfer filter only works for excitons.");
  }
  std::vector<Index> indexes;
  Lowdin low;
  QMFragment<BSE_Population> frag;

  frag.copy_withoutvalue(fragment_);
  std::vector<QMFragment<BSE_Population> > loc = {frag};
  low.CalcChargeperFragment(loc, orb, type);
  Eigen::VectorXd dq = (loc[0].value().H + loc[0].value().E).cwiseAbs();
  for (Index i = 0; i < dq.size(); i++) {
    if (dq[i] > fragment_.value()) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

void DeltaQ_filter::WriteToCpt(CheckpointWriter& w) { fragment_.WriteToCpt(w); }

void DeltaQ_filter::ReadFromCpt(CheckpointReader& r) {
  fragment_.ReadFromCpt(r);
}

}  // namespace xtp
}  // namespace votca
