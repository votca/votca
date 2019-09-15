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
#ifndef VOTCA_XTP_ESPFIT_H
#define VOTCA_XTP_ESPFIT_H

#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmfragment.h>
/**
 * \brief Takes a list of atoms, and the corresponding density matrix and puts
 * out a table of partial charges
 *
 *
 *
 */

namespace votca {
namespace xtp {
class Orbitals;
class Grid;
class QMState;
class QMMolecule;

class Espfit {
 public:
  Espfit(Logger& log) : _log(log) {
    _pairconstraint.resize(0);
    _regionconstraint.resize(0);
  }

  void setUseSVD(double conditionnumber) {
    _do_svd = true;
    _conditionnumber = conditionnumber;
  }

  void setPairConstraint(std::vector<std::pair<int, int> > pairconstraint) {
    _pairconstraint = pairconstraint;
  }

  void setRegionConstraint(std::vector<QMFragment<double> > regionconstraint) {
    _regionconstraint = regionconstraint;
  }

  StaticSegment Fit2Density(const Orbitals& orbitals, const QMState& state,
                            std::string gridsize);

 private:
  Logger& _log;
  bool _do_svd = true;
  double _conditionnumber = 1e-8;

  std::vector<std::pair<int, int> > _pairconstraint;  //  pairconstraint[i] is
                                                      //  all the atomindices
                                                      //  which have the same
                                                      //  charge

  std::vector<QMFragment<double> > _regionconstraint;

  void EvalNuclearPotential(const QMMolecule& atoms, Grid& grid);

  // Fits partial charges to Potential on a grid, constrains net charge
  StaticSegment FitPartialCharges(const Orbitals& orbitals, const Grid& grid,
                                  double netcharge);
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ESPFIT_H
