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
#ifndef VOTCA_XTP_STATICREGION_H
#define VOTCA_XTP_STATICREGION_H

#include <votca/xtp/mmregion.h>

namespace votca {
namespace xtp {
class QMRegion;
class PolarRegion;
class StaticRegion;

class StaticRegion : public MMRegion<StaticSegment> {
 public:
  StaticRegion(int id, Logger& log) : MMRegion<StaticSegment>(id, log) {}

  std::string identify() const override { return "static"; }

  void Initialize(const tools::Property& prop) override { return; }

  bool Converged() const override { return true; }

  double Etotal() const override { return 0.0; }

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions) override {
    return;
  }
  void Reset() override { return; };

 protected:
  void ResetRegion() { return; }
  void AppendResult(tools::Property& prop) const override { return; }
  double InteractwithQMRegion(const QMRegion& region) override { return 0.0; }
  double InteractwithPolarRegion(const PolarRegion& region) override {
    return 0.0;
  }
  double InteractwithStaticRegion(const StaticRegion& region) override {
    return 0.0;
  }

 private:
};

}  // namespace xtp
}  // namespace votca

#endif /* VOTCA_XTP_MMREGION_H */
