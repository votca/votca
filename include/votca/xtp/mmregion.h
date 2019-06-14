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
#ifndef VOTCA_XTP_MMREGION_H
#define VOTCA_XTP_MMREGION_H

#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/region.h>

namespace votca {
namespace xtp {

class QMRegion;
class PolarRegion;
class StaticRegion;

template <class T>
class MMRegion : public Region {
 public:
  MMRegion(int id, Logger& log) : Region(id, log){};
  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

  int size() const { return _segments.size(); }

  typedef typename std::vector<T>::iterator iterator;

  virtual void Initialize(const tools::Property& prop) = 0;

  virtual bool Converged() const = 0;

  virtual void Evaluate(std::vector<std::unique_ptr<Region> >& regions) = 0;

  virtual std::string identify() const = 0;

  const T& operator[](int index) const { return _segments[index]; }
  T& operator[](int index) { return _segments[index]; }

  typename std::vector<T>::iterator begin() { return _segments.begin(); }
  typename std::vector<T>::iterator end() { return _segments.end(); }

  typename std::vector<T>::const_iterator begin() const {
    return _segments.begin();
  }
  typename std::vector<T>::const_iterator end() const {
    return _segments.end();
  }

  void Reset();

  void WritePDB(csg::PDBWriter& writer) const;

  void push_back(const T& seg) { _segments.push_back(seg); }

 protected:
  virtual void InteractwithQMRegion(const QMRegion& region) = 0;
  virtual void InteractwithPolarRegion(const PolarRegion& region) = 0;
  virtual void InteractwithStaticRegion(const StaticRegion& region) = 0;

  std::vector<T> _segments;
};

}  // namespace xtp
}  // namespace votca

#endif /* VOTCA_XTP_MMREGION_H */
