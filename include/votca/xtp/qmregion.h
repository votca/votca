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
#include <votca/xtp/region.h>

#include "orbitals.h"

#pragma once
#ifndef VOTCA_XTP_QMREGION_H
#define VOTCA_XTP_QMREGION_H

/**
 * \brief base class to derive regions from
 *
 *
 *
 */

namespace votca {
namespace xtp {

class PolarRegion;
class StaticRegion;
class QMRegion : public Region {

 public:
  QMRegion(int id, Logger& log) : Region(id, log){};
  ~QMRegion(){};

  void Initialize(const tools::Property& prop);

  bool Converged() const;

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions);

  void ApplyInfluenceOfOtherRegions(
      const std::vector<std::unique_ptr<Region> >& regions);

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

  int size() const { return _size; }

  void WritePDB(csg::PDBWriter& writer) const;

  std::string identify() const { return "QMRegion"; }

  void push_back(const QMMolecule& mol) {
    if (_orb.QMAtoms().size() == 0) {
      _orb.QMAtoms() = mol;
    } else {
      _orb.QMAtoms().AddContainer(mol);
    }
    _size++;
  }

 protected:
  void ResetRegion();
  void InteractwithQMRegion(QMRegion& region);
  void InteractwithPolarRegion(PolarRegion& region);
  void InteractwithStaticRegion(StaticRegion& region);

 private:
  int _size = 0;
  Orbitals _orb;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_REGION_H
