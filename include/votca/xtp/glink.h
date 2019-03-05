/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_XTP_GLINK_H
#define VOTCA_XTP_GLINK_H
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {

class GNode;
class GLink {

 public:
  GLink(GNode* dest, double rate, const Eigen::Vector3d& dr)
      : destination(dest), _rate(rate), _dr(dr), _decayevent(false){};

  GLink(double rate) : _rate(rate), _decayevent(true){};

  double getValue() const { return _rate; }
  double getRate() const { return _rate; }
  GNode* getDestination() const { return destination; }
  const Eigen::Vector3d& getDeltaR() const { return _dr; }

  bool isDecayEvent() const { return _decayevent; }

 private:
  GNode* destination = nullptr;
  double _rate = 0.0;
  Eigen::Vector3d _dr = Eigen::Vector3d::Zero();
  bool _decayevent = false;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GLINK_H
