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
#ifndef VOTCA_XTP_QMSTATE_H
#define VOTCA_XTP_QMSTATE_H

// Standard includes
#include <string>

// VOTCA includes
#include <votca/tools/types.h>

namespace votca {
namespace xtp {

class QMStateType {
 public:
  enum statetype {  // never change the values
    Electron = 0,
    Hole = 1,
    Singlet = 2,
    Triplet = 3,
    Gstate = 4,
    PQPstate,
    DQPstate,
    KSstate
  };

  QMStateType(const statetype& type) : type_(type) { ; }
  QMStateType() { ; }
  QMStateType(const std::string& s) { FromString(s); }

  statetype Type() const { return type_; }

  void FromString(const std::string& statetypestring);

  std::string ToString() const;

  std::string ToLongString() const;

  bool operator==(const QMStateType& rhs) const { return type_ == rhs.Type(); }

  bool operator!=(const QMStateType& rhs) const { return type_ != rhs.Type(); }

  bool operator==(const QMStateType::statetype& rhs) const {
    return type_ == rhs;
  }

  bool operator!=(const QMStateType::statetype& rhs) const {
    return type_ != rhs;
  }

  bool isExciton() const {
    return (type_ == statetype::Singlet || type_ == statetype::Triplet);
  }

  bool isKMCState() const {
    return (type_ == statetype::Singlet || type_ == statetype::Triplet ||
            type_ == statetype::Hole || type_ == statetype::Electron);
  }

  bool isSingleParticleState() const {
    return (type_ == statetype::PQPstate || type_ == statetype::DQPstate ||
            type_ == KSstate);
  }

  bool isGWState() const {
    return (type_ == statetype::PQPstate || type_ == statetype::DQPstate);
  }

  bool isKSState() const { return (type_ == statetype::KSstate); }

  bool isPQPState() const { return (type_ == statetype::PQPstate); }

 private:
  statetype type_;
};

/**
 *  \brief  Storage class for properties of QMStateTypes, which can be used in
 * KMC
 *
 *
 */

template <class T>
class QMStateCarrierStorage {
 public:
  QMStateCarrierStorage() { content_ = {0, 0, 0, 0}; }

  void setValue(T value, QMStateType t) {
    assert(t.isKMCState() &&
           "QMStateCarrierStorage QMStateType is not for KMC simulations");
    content_[t.Type()] = value;
  }

  T getValue(QMStateType t) const {
    assert(t.isKMCState() &&
           "QMStateCarrierStorage QMStateType is not for KMC simulations");
    return content_[t.Type()];
  }

 private:
  std::array<T, 4> content_;
};

/**
 *  \brief  Identifier for QMstates. Strings like S1 are converted into enum
 * +zero indexed int
 *
 *
 */

class QMState {

 public:
  QMState(const QMStateType::statetype& type, Index index, bool transition)
      : type_(QMStateType(type)), index_(index), transition_(transition) {
    ;
  }
  QMState(const QMStateType& type, Index index, bool transition)
      : type_(type), index_(index), transition_(transition) {
    ;
  }
  QMState() { ; }
  QMState(const std::string& statestring) { FromString(statestring); }
  void FromString(const std::string& statestring);

  std::string ToString() const;

  std::string ToLongString() const;

  const QMStateType& Type() const { return type_; }

  bool isTransition() const { return transition_; }
  Index StateIdx() const { return index_; }

  bool operator==(const QMState& rhs) const {
    return (type_ == rhs.Type() && index_ == rhs.StateIdx());
  }

  bool operator!=(const QMState& rhs) const {
    return (type_ != rhs.Type() || index_ != rhs.StateIdx());
  }

 private:
  Index DetermineIndex(const std::string& statestring);
  QMStateType DetermineType(const std::string& statestring);
  QMStateType type_;

  Index index_;

  bool transition_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMSTATE_H
