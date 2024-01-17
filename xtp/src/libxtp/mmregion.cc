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

// Local VOTCA includes
#include "votca/xtp/mmregion.h"

namespace votca {
namespace xtp {

template <class T>
double MMRegion<T>::charge() const {
  double charge = 0.0;
  for (const auto& seg : segments_) {
    for (const auto& site : seg) {
      charge += site.getCharge();
    }
  }
  return charge;
}

template <class T>
void MMRegion<T>::WritePDB(csg::PDBWriter& writer) const {
  for (const auto& seg : segments_) {
    writer.WriteContainer(seg);
  }
}

template <class T>
void MMRegion<T>::WriteToCpt(CheckpointWriter& w) const {
  w(id_, "id");
  w(identify(), "type");
  Index size = Index(segments_.size());
  w(size, "size");
  CheckpointWriter ww = w.openChild("segments");
  for (const auto& seg : segments_) {
    CheckpointWriter www =
        ww.openChild(seg.identify() + "_" + std::to_string(seg.getId()));
    seg.WriteToCpt(www);
  }
}
template <class T>
void MMRegion<T>::ReadFromCpt(CheckpointReader& r) {
  r(id_, "id");
  Index size;
  r(size, "size");
  segments_.clear();
  segments_.reserve(size);
  T dummy("dummy", 0);
  CheckpointReader rr = r.openChild("segments");
  std::vector<std::string> names = rr.getChildGroupNames();
  if ( Index(names.size()) != size ){
      std::stringstream message;
      message << "Size inconsistency in region " << std::endl;
      throw std::runtime_error(message.str());
  }
  for (auto name : names ){
    CheckpointReader rrr =
        rr.openChild(name);
    segments_.push_back(T(rrr));
  }
}

template class MMRegion<PolarSegment>;
template class MMRegion<StaticSegment>;

}  // namespace xtp
}  // namespace votca
