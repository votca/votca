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
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/region.h>

#ifndef VOTCA_XTP_CLASSICALREGION_H
#define VOTCA_XTP_CLASSICALREGION_H

/**
 * \brief base class to derive regions from
 *
 *
 *
 */

namespace votca {
namespace xtp {

template <class T>
class ClassicalRegion<T> : public Region {

 public:
  ~ClassicalRegion(){};

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

  int size() const { return _segments.size(); }

  std::string identify() const {return ""};

 private:
  std::vector<T> _segments;
};
template <class T>
ClassicalRegion<T>::WriteToCpt(CheckpointWriter& w) const {
  w("name", _name);
  w("id", _id);
  w("type", identify()) for (const auto& seg : _segments) {
    w.openChild(seg.identify) seg.WriteToCpt(w);
  }
}

template <>
std::string ClassicalRegion<PolarSegment>::identify() const {
    return "PolarRegion"};
template <>
std::string ClassicalRegion<StaticRegion>::identify() const {
    return "StaticRegion"};

typedef ClassicalRegion<PolarSegment> PolarRegion;
typedef ClassicalRegion<StaticSegment> StaticRegion;

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_REGION_H
