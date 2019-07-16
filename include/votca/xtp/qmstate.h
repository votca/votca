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

#pragma once
#ifndef VOTCA_XTP_QMSTATE_H
#define VOTCA_XTP_QMSTATE_H

#include <string>

#include "checkpointreader.h"
#include "checkpointwriter.h"

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

  QMStateType(const statetype& type) : _type(type) { ; }
  QMStateType() { ; }
  QMStateType(const std::string& s) { FromString(s); }

  statetype Type() const { return _type; }

  void FromString(const std::string& statetypestring);

  std::string ToString() const;

  std::string ToLongString() const;

  bool operator==(const QMStateType& rhs) const { return _type == rhs.Type(); }

  bool operator!=(const QMStateType& rhs) const { return _type != rhs.Type(); }

  bool operator==(const QMStateType::statetype& rhs) const {
    return _type == rhs;
  }

  bool operator!=(const QMStateType::statetype& rhs) const {
    return _type != rhs;
  }

  bool isExciton() const {
    return (_type == statetype::Singlet || _type == statetype::Triplet);
  }

  bool isKMCState() const {
    return (_type == statetype::Singlet || _type == statetype::Triplet ||
            _type == statetype::Hole || _type == statetype::Electron);
  }

  bool isSingleParticleState() const {
    return (_type == statetype::PQPstate || _type == statetype::DQPstate ||
            _type == KSstate);
  }

  bool isGWState() const {
    return (_type == statetype::PQPstate || _type == statetype::DQPstate);
  }

 private:
  statetype _type;
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
  QMStateCarrierStorage() { _content = {0, 0, 0, 0}; }

  void setValue(T value, QMStateType t) {
    assert(t.isKMCState() &&
           "QMStateCarrierStorage QMStateType is not for KMC simulations");
    _content[t.Type()] = value;
  }

  T getValue(QMStateType t) const {
    assert(t.isKMCState() &&
           "QMStateCarrierStorage QMStateType is not for KMC simulations");
    return _content[t.Type()];
  }

 private:
  std::array<T, 4> _content;
};

/**
 *  \brief  Identifier for QMstates. Strings like S1 are converted into enum
 * +zero indexed int
 *
 *
 */

class QMState {

 public:
  QMState(const QMStateType::statetype& type, int index, bool transition)
      : _type(QMStateType(type)), _index(index), _transition(transition) {
    ;
  }
  QMState(const QMStateType& type, int index, bool transition)
      : _type(type), _index(index), _transition(transition) {
    ;
  }
  QMState() { ; }
  QMState(const std::string& statestring) { FromString(statestring); }
  void FromString(const std::string& statestring);

  std::string ToString() const;

  std::string ToLongString() const;

  const QMStateType& Type() const { return _type; }

  bool isTransition() const { return _transition; }
  int Index() const { return _index; }

  bool operator==(const QMState& rhs) const {
    return (_type == rhs.Type() && _index == rhs.Index());
  }

  bool operator!=(const QMState& rhs) const {
    return (_type != rhs.Type() || _index != rhs.Index());
  }

 private:
  int DetermineIndex(const std::string& statestring);
  QMStateType DetermineType(const std::string& statestring);
  QMStateType _type;

  int _index;

  bool _transition;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMSTATE_H
