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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#pragma once
#ifndef VOTCA_XTP_QMPAIR_H
#define VOTCA_XTP_QMPAIR_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "qmstate.h"
#include "segment.h"

namespace votca {
namespace xtp {

class QMPair {
 public:
  enum PairType { Hopping = 0, Excitoncl = 1 };

  static std::string get_name(PairType type) {
    switch (type) {
      case Hopping:
        return "Hopping";
      case Excitoncl:
        return "Excitoncl";
        // no default case to trigger compiler error
    }
    return "";
  }

  struct data {
    Index id;
    Index Seg1Id;
    Index Seg2Id;
    double RX;
    double RY;
    double RZ;

    char* pair_type;

    double lambda0e;
    double lambda0h;
    double lambda0s;
    double lambda0t;

    double jeff2e;
    double jeff2h;
    double jeff2s;
    double jeff2t;
  };

  static PairType get_Enum(std::string type) {
    if (type == "Hopping") {
      return PairType::Hopping;
    } else if (type == "Excitoncl") {
      return PairType::Excitoncl;
    } else {
      throw std::runtime_error("get_Enum input is invalid");
    }
  }

  QMPair(Index id, const Segment* seg1, const Segment* seg2,
         const Eigen::Vector3d& delta_R);

  QMPair(const data& d, const std::vector<Segment>& segments) {
    ReadData(d, segments);
  }

  Index getId() const { return id_; }
  void setId(Index id) { id_ = id; }

  const Eigen::Vector3d& R() const { return R_; }
  double Dist() const { return R_.norm(); }

  void setLambdaO(double lO, QMStateType state) {
    lambda0_.setValue(lO, state);
  }
  double getLambdaO(QMStateType state) const {
    return lambda0_.getValue(state);
  }

  double getReorg12(QMStateType state) const {
    return segments_.first->getU_nX_nN(state) +
           segments_.second->getU_xN_xX(state);
  }  // 1->2
  double getReorg21(QMStateType state) const {
    return segments_.first->getU_xN_xX(state) +
           segments_.second->getU_nX_nN(state);
  }  // 2->1

  double getJeff2(QMStateType state) const { return Jeff2_.getValue(state); }
  void setJeff2(double Jeff2, QMStateType state) {
    Jeff2_.setValue(Jeff2, state);
  }

  double getdE12(QMStateType state) const {
    return segments_.first->getSiteEnergy(state) -
           segments_.second->getSiteEnergy(state);
  }

  Segment Seg2PbCopy() const;
  const Segment* Seg1() const { return segments_.first; }
  const Segment* Seg2() const { return segments_.second; }

  const Segment* first() { return segments_.first; }
  const Segment* second() { return segments_.second; }

  void setType(PairType pair_type) { pair_type_ = pair_type; }
  const PairType& getType() const { return pair_type_; }

  static void SetupCptTable(CptTable& table);
  void WriteData(data& d) const;

  void ReadData(const data& d, const std::vector<Segment>& segments);

 private:
  Index id_ = -1;
  std::pair<const Segment*, const Segment*> segments_;

  Eigen::Vector3d R_ = Eigen::Vector3d::Zero();

  PairType pair_type_ = PairType::Hopping;

  QMStateCarrierStorage<double> lambda0_;
  QMStateCarrierStorage<double> Jeff2_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMPAIR_H
