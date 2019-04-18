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

#include <votca/xtp/mmregion.h>

namespace votca {
namespace xtp {
template <class T>
std::string MMRegion<T>::identify() const {
  return "";
}

template <>
std::string PolarRegion::identify() const {
  return "PolarRegion";
}

template <>
std::string StaticRegion::identify() const {
  return "StaticRegion";
}

template <class T>
void MMRegion<T>::WriteToCpt(CheckpointWriter& w) const {
  w(_name, "name");
  w(_id, "id");
  w(identify(), "type");
  for (const auto& seg : _segments) {
    w.openChild(seg.identify() + "_" + std::to_string(seg.getId()));
    seg.WriteToCpt(w);
  }
}
template <class T>
void MMRegion<T>::ReadFromCpt(CheckpointReader& r) {}

template class MMRegion<PolarSegment>;
template class MMRegion<StaticSegment>;

}  // namespace xtp
}  // namespace votca
