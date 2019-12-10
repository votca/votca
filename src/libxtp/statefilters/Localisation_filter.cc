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

#include "Localisation_filter.h"
#include <votca/xtp/populationanalysis.h>
namespace votca {
namespace xtp {

void Localisation_filter::Initialize(const tools::Property& options) {
  std::string indices =
      options.ifExistsReturnElseThrowRuntimeError<std::string>("fragment");
  _fragment = QMFragment<double>(0, indices);
  _fragment.value() =
      options.ifExistsReturnElseThrowRuntimeError<double>("threshold");
}

void Localisation_filter::Info(Logger& log) const {
  XTP_LOG(Log::error, log) << "Using localisation tracker for fragment "
                           << _fragment << std::flush;
}

void Localisation_filter::UpdateHist(const Orbitals&, QMState) { return; }

std::vector<Index> Localisation_filter::CalcIndeces(const Orbitals& orb,
                                                    QMStateType type) const {

  if (!type.isExciton()) {
    throw std::runtime_error("Localisation filter only works for excitons.");
  }
  std::vector<Index> indexes;
  Lowdin low;
  QMFragment<BSE_Population> frag;
  frag.copy_withoutvalue(_fragment);
  std::vector<QMFragment<BSE_Population> > loc = {frag};
  low.CalcChargeperFragment(loc, orb, type);
  const Eigen::VectorXd& popE = loc[0].value().E;
  const Eigen::VectorXd& popH = loc[0].value().H;
  for (Index i = 0; i < popE.size(); i++) {
    if (std::abs(popE[i]) > _fragment.value() &&
        std::abs(popH[i]) > _fragment.value()) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

void Localisation_filter::WriteToCpt(CheckpointWriter& w) {
  _fragment.WriteToCpt(w);
}

void Localisation_filter::ReadFromCpt(CheckpointReader& r) {
  _fragment.ReadFromCpt(r);
}

}  // namespace xtp
}  // namespace votca