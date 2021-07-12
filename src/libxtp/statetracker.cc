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
#include "votca/xtp/statetracker.h"
#include "votca/xtp/filterfactory.h"

namespace votca {
namespace xtp {
using std::flush;

void StateTracker::Initialize(const tools::Property& options) {

  FilterFactory::RegisterAll();
  for (const tools::Property& filter : options) {
    filters_.push_back(Filter().Create(filter.name()));
  }

  for (auto& filter : filters_) {
    const tools::Property& filterop = options.get(filter->Identify());
    filter->Initialize(filterop);
  }
}

void StateTracker::PrintInfo() const {
  XTP_LOG(Log::error, *log_)
      << "Initial state: " << statehist_[0].ToString() << flush;
  if (statehist_.size() > 1) {
    XTP_LOG(Log::error, *log_)
        << "Last state: " << statehist_.back().ToString() << flush;
  }

  if (filters_.empty()) {
    XTP_LOG(Log::error, *log_) << "WARNING: No tracker is used " << flush;
  } else {
    for (const auto& filter : filters_) {
      filter->Info(*log_);
    }
  }
}

std::vector<Index> StateTracker::ComparePairofVectors(
    std::vector<Index>& vec1, std::vector<Index>& vec2) const {
  std::vector<Index> result(std::min(vec1, vec2));
  std::sort(vec1.begin(), vec1.end());
  std::sort(vec2.begin(), vec2.end());
  std::vector<Index>::iterator it = std::set_intersection(
      vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), result.begin());
  result.resize(it - result.begin());
  return result;
}

std::vector<Index> StateTracker::CollapseResults(
    std::vector<std::vector<Index>>& results) const {
  if (results.empty()) {
    return std::vector<Index>(0);
  } else {
    std::vector<Index> result = results[0];
    for (Index i = 1; i < Index(results.size()); i++) {
      result = ComparePairofVectors(result, results[i]);
    }
    return result;
  }
}

QMState StateTracker::CalcState(const Orbitals& orbitals) const {

  if (filters_.empty()) {
    return statehist_[0];
  }

  std::vector<std::vector<Index>> results;
  for (const auto& filter : filters_) {
    if (statehist_.size() < 2 && filter->NeedsInitialState()) {
      XTP_LOG(Log::error, *log_)
          << "Filter " << filter->Identify()
          << " not used in first iteration as it needs a reference state"
          << flush;
      continue;
    }
    results.push_back(filter->CalcIndeces(orbitals, statehist_[0].Type()));
  }

  std::vector<Index> result = CollapseResults(results);
  QMState state;
  if (result.size() < 1) {
    state = statehist_.back();
    XTP_LOG(Log::error, *log_)
        << "No State found by tracker using last state: " << state.ToString()
        << flush;
  } else {
    state = QMState(statehist_.back().Type(), result[0], false);
    XTP_LOG(Log::error, *log_)
        << "Next State is: " << state.ToString() << flush;
  }
  return state;
}

QMState StateTracker::CalcStateAndUpdate(const Orbitals& orbitals) {
  QMState result = CalcState(orbitals);
  statehist_.push_back(result);
  for (auto& filter : filters_) {
    filter->UpdateHist(orbitals, result);
  }
  return result;
}

void StateTracker::WriteToCpt(CheckpointWriter& w) const {
  std::vector<std::string> statehiststring;
  statehiststring.reserve(statehist_.size());
  for (const QMState& s : statehist_) {
    statehiststring.push_back(s.ToString());
  }
  w(statehiststring, "statehist");

  for (const auto& filter : filters_) {
    CheckpointWriter ww = w.openChild(filter->Identify());
    filter->WriteToCpt(ww);
  }
}

void StateTracker::ReadFromCpt(CheckpointReader& r) {
  FilterFactory::RegisterAll();
  std::vector<std::string> statehiststring;
  r(statehiststring, "statehist");
  statehist_.clear();
  statehist_.reserve(statehiststring.size());
  for (const std::string& s : statehiststring) {
    statehist_.push_back(QMState(s));
  }
  filters_.clear();
  for (const std::string filtername : r.getChildGroupNames()) {

    CheckpointReader rr = r.openChild(filtername);
    filters_.push_back(Filter().Create(filtername));
    filters_.back()->ReadFromCpt(rr);
  }
}

}  // namespace xtp
}  // namespace votca
