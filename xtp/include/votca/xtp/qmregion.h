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

#pragma once
#ifndef VOTCA_XTP_QMREGION_H
#define VOTCA_XTP_QMREGION_H

// Local VOTCA includes
#include "hist.h"
#include "orbitals.h"
#include "qmpackagefactory.h"
#include "region.h"
#include "statetracker.h"
#include "logger.h"

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
  QMRegion(Index id, Logger& log, std::string workdir)
      : Region(id, log), workdir_(workdir) {
    QMPackageFactory::RegisterAll();
  };
  ~QMRegion() override = default;

  void Initialize(const tools::Property& prop) override;

  bool Converged() const override;

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions) override;

  void WriteToCpt(CheckpointWriter& w) const override;

  void ReadFromCpt(CheckpointReader& r) override;

  void ApplyQMFieldToPolarSegments(std::vector<PolarSegment>& segments) const;

  Index size() const override { return size_; }

  void WritePDB(csg::PDBWriter& writer) const override;

  std::string identify() const override { return "qmregion"; }

  void push_back(const QMMolecule& mol);

  void Reset() override;

  double charge() const override;
  double Etotal() const override { return E_hist_.back(); }

 protected:
  void AppendResult(tools::Property& prop) const override;
  double InteractwithQMRegion(const QMRegion& region) override;
  double InteractwithPolarRegion(const PolarRegion& region) override;
  double InteractwithStaticRegion(const StaticRegion& region) override;

 private:
  void AddNucleiFields(std::vector<PolarSegment>& segments,
                       const StaticSegment& seg) const;

  Index size_ = 0;
  Orbitals orb_;

  QMState initstate_;
  std::string workdir_ = "";
  std::unique_ptr<QMPackage> qmpackage_ = nullptr;

  std::string grid_accuracy_for_ext_interaction_ = "medium";

  hist<double> E_hist_;
  hist<Eigen::MatrixXd> Dmat_hist_;

  // convergence options
  double DeltaD_ = 5e-5;
  double DeltaE_ = 5e-5;
  double DeltaDmax_ = 5e-5;

  bool do_gwbse_ = false;
  bool do_localize_ = false;
  bool do_dft_in_dft_ = false;

  tools::Property dftoptions_;
  tools::Property gwbseoptions_;
  tools::Property localize_options_;

  StateTracker statetracker_;

  //Logger* logger_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMREGION_H
