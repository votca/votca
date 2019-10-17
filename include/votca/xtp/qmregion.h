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
#ifndef VOTCA_XTP_QMREGION_H
#define VOTCA_XTP_QMREGION_H

#include <votca/xtp/region.h>

#include "orbitals.h"
#include "statetracker.h"
#include <votca/xtp/hist.h>
#include <votca/xtp/qmpackagefactory.h>
/**
 * \brief defines a qm region and runs dft and gwbse calculations
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
  QMRegion(int id, Logger& log, std::string workdir)
      : Region(id, log), _workdir(workdir) {
    QMPackageFactory::RegisterAll();
  };
  ~QMRegion() override = default;
  ;

  void Initialize(const tools::Property& prop) override;

  bool Converged() const override;

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions) override;

  void WriteToCpt(CheckpointWriter& w) const override;

  void ReadFromCpt(CheckpointReader& r) override;

  void ApplyQMFieldToPolarSegments(std::vector<PolarSegment>& segments) const;

  int size() const override { return _size; }

  void WritePDB(csg::PDBWriter& writer) const override;

  std::string identify() const override { return "qm"; }

  void push_back(const QMMolecule& mol);

  void Reset() override;

  double charge() const override;
  double Etotal() const override { return _E_hist.back(); }

 protected:
  void AppendResult(tools::Property& prop) const override;
  double InteractwithQMRegion(const QMRegion& region) override;
  double InteractwithPolarRegion(const PolarRegion& region) override;
  double InteractwithStaticRegion(const StaticRegion& region) override;

 private:
  void AddNucleiFields(std::vector<PolarSegment>& segments,
                       const StaticSegment& seg) const;

  int _size = 0;
  Orbitals _orb;

  QMState _initstate;
  std::string _workdir = "";
  std::unique_ptr<QMPackage> _qmpackage = nullptr;

  std::string _grid_accuracy_for_ext_interaction = "medium";

  hist<double> _E_hist;
  hist<Eigen::MatrixXd> _Dmat_hist;

  // convergence options
  double _DeltaD = 5e-5;
  double _DeltaE = 5e-5;

  bool _do_gwbse = false;

  tools::Property _dftoptions;
  tools::Property _gwbseoptions;

  StateTracker _statetracker;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_REGION_H
