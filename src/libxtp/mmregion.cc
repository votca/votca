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
void MMRegion<T>::WritePDB(csg::PDBWriter& writer) const {
  for (const auto& seg : _segments) {
    writer.WriteContainer(seg);
  }
}

template <class T>
void MMRegion<T>::Reset() {
  XTP_LOG_SAVE(logINFO, _log)
      << TimeStamp() << " Removed all previous values from region"
      << std::flush;
  for (auto& seg : _segments) {
    for (auto& site : seg) {
      site.Reset();
    }
  }
}

template <class T>
void MMRegion<T>::WriteToCpt(CheckpointWriter& w) const {
  w(_id, "id");
  w(identify(), "type");
  int size = _segments.size();
  w(size, "size");
  CheckpointWriter ww = w.openChild("segments");
  for (const auto& seg : _segments) {
    CheckpointWriter www =
        ww.openChild(seg.identify() + "_" + std::to_string(seg.getId()));
    seg.WriteToCpt(www);
  }
}
template <class T>
void MMRegion<T>::ReadFromCpt(CheckpointReader& r) {
  r(_id, "id");
  int size;
  r(size, "size");
  _segments.clear();
  _segments.reserve(size);
  T dummy("dummy", 0);
  CheckpointReader rr = r.openChild("segments");
  for (int i = 0; i < size; i++) {
    CheckpointReader rrr =
        rr.openChild(dummy.identify() + "_" + std::to_string(i));
    _segments.push_back(T(rrr));
  }
}

template class MMRegion<PolarSegment>;
template class MMRegion<StaticSegment>;

}  // namespace xtp
}  // namespace votca
