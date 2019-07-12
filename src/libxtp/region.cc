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

std::vector<double> Region::ApplyInfluenceOfOtherRegions(
    std::vector<std::unique_ptr<Region> >& regions) {
  std::vector<double> energies = std::vector<double>(regions.size(), 0.0);
  for (std::unique_ptr<Region>& reg : regions) {
    int id = reg->getId();
    if (id == this->getId()) {
      continue;
    }
    if (!reg->Evaluated()) {
      XTP_LOG_SAVE(logINFO, _log)
          << TimeStamp() << " Skipping interaction between " << this->identify()
          << " " << this->getId() << " and " << reg->identify() << " "
          << reg->getId() << std::flush;
      continue;
    }

    QMRegion QMdummy(0, _log, "");
    StaticRegion Staticdummy(0, _log);
    PolarRegion Polardummy(0, _log);
    XTP_LOG_SAVE(logINFO, _log)
        << TimeStamp() << " Evaluating interaction between " << this->identify()
        << " " << this->getId() << " and " << reg->identify() << " "
        << reg->getId() << std::flush;
    if (reg->identify() == QMdummy.identify()) {
      QMRegion* qmregion = dynamic_cast<QMRegion*>(reg.get());
      energies[id] = InteractwithQMRegion(*qmregion);
    } else if (reg->identify() == Staticdummy.identify()) {
      StaticRegion* staticregion = dynamic_cast<StaticRegion*>(reg.get());
      energies[id] = InteractwithStaticRegion(*staticregion);
    } else if (reg->identify() == Polardummy.identify()) {
      PolarRegion* polarregion = dynamic_cast<PolarRegion*>(reg.get());
      energies[id] = InteractwithPolarRegion(*polarregion);
    } else {
      throw std::runtime_error(
          "Interaction of regions with types:" + this->identify() + " and " +
          reg->identify() + " not implemented");
    }
  }

  return energies;
}

void Region::AddResults(tools::Property& prop) const {
  tools::Property& region = prop.add("region", "");
  region.setAttribute("type", identify());
  region.setAttribute("id", getId());
  region.setAttribute("size", size());
  region.setAttribute("Tot_charge", (boost::format("%1$1.6e") % charge()));
  AppendResult(region);
}

}  // namespace xtp
}  // namespace votca
