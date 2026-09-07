/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

// Standard includes
#include <sstream>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/ewaldregistry.h"

namespace votca {
namespace xtp {

namespace {

std::string StateTag(EwaldChargeState state) {
  switch (state) {
    case EwaldChargeState::Neutral:
      return "neutral";
    case EwaldChargeState::Electron:
      return "electron";
    case EwaldChargeState::Hole:
      return "hole";
  }
  // Unreachable for a valid enum value; kept explicit rather than relying
  // on undefined behaviour if the enum is ever extended without updating
  // this function.
  throw std::runtime_error("EwaldRegistry: unknown EwaldChargeState");
}

bool TagToState(const std::string& tag, EwaldChargeState& state) {
  if (tag == "neutral") {
    state = EwaldChargeState::Neutral;
    return true;
  }
  if (tag == "electron") {
    state = EwaldChargeState::Electron;
    return true;
  }
  if (tag == "hole") {
    state = EwaldChargeState::Hole;
    return true;
  }
  return false;
}

std::string GroupName(Index id, EwaldChargeState state) {
  return std::to_string(id) + "_" + StateTag(state);
}

}  // namespace

void EwaldRegistry::Register(Index id, EwaldChargeState state,
                              PolarSegment segment) {
  // PolarSegment has no default constructor, so operator[] (which would
  // value-initialize a fresh entry before assigning into it) is not usable
  // here; insert_or_assign constructs the entry directly from the moved
  // argument instead.
  store_.insert_or_assign(Key(id, state), std::move(segment));
}

bool EwaldRegistry::Has(Index id, EwaldChargeState state) const {
  return store_.find(Key(id, state)) != store_.end();
}

const PolarSegment& EwaldRegistry::Get(Index id, EwaldChargeState state) const {
  auto it = store_.find(Key(id, state));
  if (it == store_.end()) {
    std::stringstream message;
    message << "EwaldRegistry: no entry for segment id " << id
            << " in state '" << StateTag(state) << "'";
    throw std::out_of_range(message.str());
  }
  return it->second;
}

PolarSegment& EwaldRegistry::Get(Index id, EwaldChargeState state) {
  auto it = store_.find(Key(id, state));
  if (it == store_.end()) {
    std::stringstream message;
    message << "EwaldRegistry: no entry for segment id " << id
            << " in state '" << StateTag(state) << "'";
    throw std::out_of_range(message.str());
  }
  return it->second;
}

void EwaldRegistry::Erase(Index id, EwaldChargeState state) {
  store_.erase(Key(id, state));
}

std::vector<Index> EwaldRegistry::AllIds() const {
  std::vector<Index> ids;
  for (const auto& entry : store_) {
    Index id = entry.first.first;
    if (ids.empty() || ids.back() != id) {
      ids.push_back(id);
    }
  }
  return ids;
}

std::vector<EwaldChargeState> EwaldRegistry::StatesFor(Index id) const {
  std::vector<EwaldChargeState> states;
  for (const auto& entry : store_) {
    if (entry.first.first == id) {
      states.push_back(entry.first.second);
    }
  }
  return states;
}

void EwaldRegistry::WriteToCpt(CheckpointWriter& w) const {
  Index size = static_cast<Index>(store_.size());
  w(size, "size");
  CheckpointWriter ww = w.openChild("entries");
  for (const auto& entry : store_) {
    Index id = entry.first.first;
    EwaldChargeState state = entry.first.second;
    const PolarSegment& seg = entry.second;
    CheckpointWriter www = ww.openChild(GroupName(id, state));
    seg.WriteToCpt(www);
  }
}

void EwaldRegistry::ReadFromCpt(CheckpointReader& r) {
  Index size;
  r(size, "size");
  store_.clear();
  CheckpointReader rr = r.openChild("entries");
  std::vector<std::string> names = rr.getChildGroupNames();
  if (Index(names.size()) != size) {
    std::stringstream message;
    message << "EwaldRegistry: size inconsistency reading checkpoint ("
            << names.size() << " groups found, " << size << " expected)";
    throw std::runtime_error(message.str());
  }
  for (const std::string& name : names) {
    // Group names are "<id>_<statetag>"; split on the last underscore
    // rather than the first, since the id portion never contains one.
    std::size_t split = name.rfind('_');
    if (split == std::string::npos) {
      throw std::runtime_error(
          "EwaldRegistry: malformed checkpoint group name '" + name + "'");
    }
    Index id = std::stoi(name.substr(0, split));
    std::string tag = name.substr(split + 1);
    EwaldChargeState state;
    if (!TagToState(tag, state)) {
      throw std::runtime_error(
          "EwaldRegistry: unknown charge state tag '" + tag +
          "' in checkpoint group name '" + name + "'");
    }
    CheckpointReader rrr = rr.openChild(name);
    store_.insert_or_assign(Key(id, state), PolarSegment(rrr));
  }
}

}  // namespace xtp
}  // namespace votca
