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

#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/region.h>
#include <votca/xtp/staticregion.h>

namespace votca {
namespace xtp {

void Region::ApplyInfluenceOfOtherRegions(
    std::vector<std::unique_ptr<Region> >& regions) {
  for (std::unique_ptr<Region>& reg : regions) {
    if (reg->getId() == this->getId()) {
      continue;
    }

    QMRegion QMdummy(0, _log);
    StaticRegion Staticdummy(0, _log);
    PolarRegion Polardummy(0, _log);
    XTP_LOG_SAVE(logINFO, _log)
        << "Evaluating interaction between:" << this->identify() << " "
        << this->getId() << " and " << reg->identify() << " " << reg->getId()
        << std::flush;
    if (reg->identify() == QMdummy.identify()) {
      QMRegion* qmregion = dynamic_cast<QMRegion*>(reg.get());
      InteractwithQMRegion(*qmregion);
    } else if (reg->identify() == Staticdummy.identify()) {
      StaticRegion* staticregion = dynamic_cast<StaticRegion*>(reg.get());
      InteractwithStaticRegion(*staticregion);
    } else if (reg->identify() == Polardummy.identify()) {
      PolarRegion* polarregion = dynamic_cast<PolarRegion*>(reg.get());
      InteractwithPolarRegion(*polarregion);
    } else {
      throw std::runtime_error(
          "Interaction of regions with types:" + this->identify() + " and " +
          reg->identify() + " not implemented");
    }
  }
  return;
}

}  // namespace xtp
}  // namespace votca
