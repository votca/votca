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
#ifndef VOTCA_XTP_CLASSICALSEGMENT_H
#define VOTCA_XTP_CLASSICALSEGMENT_H

// Local VOTCA includes
#include "atomcontainer.h"
#include "chargeresponse.h"
#include "polarsite.h"
#include "staticsite.h"

namespace votca {
namespace xtp {
template <class T>
class ClassicalSegment : public AtomContainer<T> {
 public:
  ClassicalSegment(std::string name, Index id) : AtomContainer<T>(name, id) {};

  ClassicalSegment(CheckpointReader& r) : AtomContainer<T>(r) {
    ReadChargeResponseFromCpt(r);
  };

  void LoadFromFile(std::string filename);

  void WriteMPS(std::string filename, std::string header) const;

  double CalcTotalQ() const;

  std::string identify() const;

  Eigen::Vector3d CalcDipole() const;

  friend std::ostream& operator<<(std::ostream& out,
                                  const ClassicalSegment<T>& container) {
    out << container.getId() << " " << container.identify() << " "
        << container.getType() << "\n";
    for (const T& atom : container) {
      out << atom;
    }
    out << std::endl;
    return out;
  }

  // Charge-response model (polar segments only): induced atomic charges
  // dq = -Pi phi next to the induced dipoles of the sites. Without it the
  // segment is an ordinary Thole segment.
  bool hasChargeResponse() const { return !chargeresponse_.empty(); }
  const ChargeResponse& getChargeResponse() const { return chargeresponse_; }
  // Pi must be in the atom order of this segment (checked by element)
  void setChargeResponse(const ChargeResponse& cr);

  void WriteToCpt(CheckpointWriter& w) const override;
  void ReadFromCpt(CheckpointReader& r) override;

  // additional access functions and tools for Ewald compatibility
  bool IsCharged() const;
  bool IsPolarizable() const;
  Eigen::Vector3d CalcStaticDipole() const;

 private:
  void ReadChargeResponseFromCpt(CheckpointReader& r);
  ChargeResponse chargeresponse_;
};

using PolarSegment = ClassicalSegment<PolarSite>;
using StaticSegment = ClassicalSegment<StaticSite>;

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CLASSICALSEGMENT_H
