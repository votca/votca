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

#include "votca/xtp/statetracker.h"
#include <votca/xtp/filterfactory.h>

namespace votca {
namespace xtp {
using std::flush;

void StateTracker::Initialize(const tools::Property& options) {

  std::string filters =
      options.ifExistsReturnElseThrowRuntimeError<std::string>("filters");
  tools::Tokenizer tok(filters, " ,;\n");
  std::vector<std::string> list_filters = tok.ToVector();

  FilterFactory::RegisterAll();
  for (const std::string filtername : list_filters) {
    _filters.push_back(
        std::unique_ptr<StateFilter_base>(Filter().Create(filtername)));
  }

  for (auto& filter : _filters) {
    const tools::Property& filterop = options.get(filter->Identify());
    filter->Initialize(filterop);
  }
}

void StateTracker::PrintInfo() const {
  XTP_LOG(Log::error, *_log)
      << "Initial state: " << _statehist[0].ToString() << flush;
  if (_statehist.size() > 1) {
    XTP_LOG(Log::error, *_log)
        << "Last state: " << _statehist.back().ToString() << flush;
  }

  if (_filters.empty()) {
    XTP_LOG(Log::error, *_log) << "WARNING: No tracker is used " << flush;
  } else {
    for (const auto& filter : _filters) {
      filter->Info(*_log);
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
    std::vector<std::vector<Index> >& results) const {
  if (results.size() == 1) {
    return results[0];
  } else {
    std::vector<Index> result = results[0];
    for (Index i = 1; i < Index(results.size()); i++) {
      result = ComparePairofVectors(result, results[i]);
    }
    return result;
  }
}

QMState StateTracker::CalcState(const Orbitals& orbitals) const {

  if (_filters.empty()) {
    return _statehist[0];
  }

  std::vector<std::vector<Index> > results;
  for (const auto& filter : _filters) {
    if (_statehist.size() < 2 && filter->NeedsInitialState()) {
      XTP_LOG(Log::error, *_log)
          << "Filter " << filter->Identify()
          << " not used in first iteration as it needs a reference state"
          << flush;
      continue;
    }
    results.push_back(filter->CalcIndeces(orbitals, _statehist[0].Type()));
  }

  std::vector<Index> result = CollapseResults(results);
  QMState state;
  if (result.size() < 1) {
    state = _statehist.back();
    XTP_LOG(Log::error, *_log)
        << "No State found by tracker using last state: " << state.ToString()
        << flush;
  } else {
    state = QMState(_statehist.back().Type(), result[0], false);
    XTP_LOG(Log::error, *_log)
        << "Next State is: " << state.ToString() << flush;
  }
  return state;
}

QMState StateTracker::CalcStateAndUpdate(const Orbitals& orbitals) {
  QMState result = CalcState(orbitals);
  for (auto& filter : _filters) {
    filter->UpdateHist(orbitals, result);
  }
  return result;
}

void StateTracker::WriteToCpt(CheckpointWriter& w) const {
  std::vector<std::string> statehiststring;
  statehiststring.reserve(_statehist.size());
  for (const QMState& s : _statehist) {
    statehiststring.push_back(s.ToString());
  }
  w(statehiststring, "statehist");

  for (const auto& filter : _filters) {
    CheckpointWriter ww = w.openChild(filter->Identify());
    filter->WriteToCpt(ww);
  }
}

void StateTracker::ReadFromCpt(CheckpointReader& r) {
  std::vector<std::string> statehiststring;
  r(statehiststring, "statehist");
  _statehist.clear();
  _statehist.reserve(statehiststring.size());
  for (const std::string& s : statehiststring) {
    _statehist.push_back(QMState(s));
  }
  for (auto& filter : _filters) {
    CheckpointReader rr = r.openChild(filter->Identify());
    filter->ReadFromCpt(rr);
  }
}

}  // namespace xtp
}  // namespace votca
