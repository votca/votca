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

// Standard includes
#include <regex>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

// Local VOTCA includes
#include "votca/xtp/qmstate.h"

namespace votca {
namespace xtp {

std::string QMStateType::ToString() const {
  std::string identifier = "";
  switch (type_) {
    case QMStateType::Singlet:
      identifier = "s";
      break;
    case QMStateType::Triplet:
      identifier = "t";
      break;
    case QMStateType::PQPstate:
      identifier = "pqp";
      break;
    case QMStateType::DQPstate:
      identifier = "dqp";
      break;
    case QMStateType::KSstate:
      identifier = "ks";
      break;
    case QMStateType::Gstate:
      identifier = "n";
      break;
    case QMStateType::Hole:
      identifier = "h";
      break;
    case QMStateType::Electron:
      identifier = "e";
      break;
    case QMStateType::LMOstate:
      identifier = "l";
      break;
  }
  return identifier;
}

std::string QMStateType::ToLongString() const {
  std::string identifier = "";
  switch (type_) {
    case QMStateType::Singlet:
      identifier = "singlet";
      break;
    case QMStateType::Triplet:
      identifier = "triplet";
      break;
    case QMStateType::PQPstate:
      identifier = "pert-qparticle";
      break;
    case QMStateType::DQPstate:
      identifier = "diag-qparticle";
      break;
    case QMStateType::KSstate:
      identifier = "Kohn-Sham-orbital";
      break;
    case QMStateType::Gstate:
      identifier = "groundstate";
      break;
    case QMStateType::Hole:
      identifier = "hole";
      break;
    case QMStateType::Electron:
      identifier = "electron";
      break;
    case QMStateType::LMOstate:
      identifier = "localized-orbital";
      break;
  }
  return identifier;
}

void QMStateType::FromString(const std::string& statetypestring) {
  std::string lower = boost::algorithm::to_lower_copy(statetypestring);
  boost::trim(lower);
  if (lower == "s" || lower == "singlet") {
    type_ = QMStateType::Singlet;
  } else if (lower == "t" || lower == "triplet") {
    type_ = QMStateType::Triplet;
  } else if (lower == "pqp" || lower == "pert-qparticle") {
    type_ = QMStateType::PQPstate;
  } else if (lower == "dqp" || lower == "diag-qparticle" || lower == "qpdiag") {
    type_ = QMStateType::DQPstate;
  } else if (lower == "ks" || lower == "kohn-sham-orbital") {
    type_ = QMStateType::KSstate;
  } else if (lower == "n" || lower == "groundstate" || lower == "gs") {
    type_ = QMStateType::Gstate;
  } else if (lower == "h" || lower == "hole") {
    type_ = QMStateType::Hole;
  } else if (lower == "e" || lower == "electron") {
    type_ = QMStateType::Electron;
  } else if (lower == "l" || lower == "localized-orbital") {
    type_ = QMStateType::LMOstate;
  } else {
    throw std::runtime_error("Statetype:" + statetypestring +
                             " not recognized");
  }
}

std::string QMState::ToLongString() const {
  Index index = index_;
  if (type_ == QMStateType::Singlet || type_ == QMStateType::Triplet) {
    index++;
  } else if (type_ == QMStateType::Gstate || type_ == QMStateType::Electron ||
             type_ == QMStateType::Hole) {
    return type_.ToLongString();
  }
  std::string result =
      type_.ToLongString() + (boost::format(" %i") % index).str();
  if (transition_) {
    result = "Groundstate to " + result;
  }
  return result;
}

std::string QMState::ToString() const {
  Index index = index_;
  if (type_ == QMStateType::Singlet || type_ == QMStateType::Triplet) {
    index++;
  } else if (type_ == QMStateType::Gstate || type_ == QMStateType::Electron ||
             type_ == QMStateType::Hole) {
    return type_.ToString();
  }
  std::string result = type_.ToString() + (boost::format("%i") % index).str();
  if (transition_) {
    result = "n2" + result;
  }
  return result;
}

Index QMState::DetermineIndex(const std::string& statestring) {

  std::smatch search;
  std::regex reg("[0-9]+");

  bool found_integer = std::regex_search(statestring, search, reg);
  if (!found_integer) {

    if (type_ == QMStateType::Hole || type_ == QMStateType::Electron) {
      return 0;
    }

    throw std::runtime_error("Found no index in string: " + statestring);
  }
  if (search.size() > 1) {
    throw std::runtime_error("Found more than 1 index in string: " +
                             statestring);
  }

  Index index = boost::lexical_cast<Index>(search.str(0));
  if (type_.isExciton() || type_ == QMStateType::Electron ||
      type_ == QMStateType::Hole) {
    index--;
  }
  return index;
}

QMStateType QMState::DetermineType(const std::string& statestring) {
  std::regex reg("[^0-9]+");
  std::smatch search;

  bool found_typestring = std::regex_search(statestring, search, reg);
  if (!found_typestring) {
    throw std::runtime_error("Found no type in string: " + statestring);
  }
  if (search.size() > 1) {
    throw std::runtime_error("Found more than one type in string: " +
                             statestring);
  }
  return QMStateType(search.str(0));
}

void QMState::FromString(const std::string& statestring) {
  std::string lower = boost::algorithm::to_lower_copy(statestring);
  boost::trim(lower);
  std::string rest;
  if (boost::starts_with(lower, "n2")) {
    transition_ = true;
    rest = lower.substr(2);
  } else if (boost::starts_with(lower, "groundstate to")) {
    transition_ = true;
    rest = lower.substr(14);
  } else {
    rest = lower;
    transition_ = false;
  }
  boost::trim(rest);

  type_ = DetermineType(rest);
  if (type_ != QMStateType::Singlet && transition_ == true) {
    throw std::runtime_error("Transition states only exist for singlets.");
  }
  if (type_ == QMStateType::Gstate) {
    index_ = -1;
  } else {
    index_ = DetermineIndex(rest);
  }
}

}  // namespace xtp
}  // namespace votca
