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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_SEGMENT_H
#define VOTCA_XTP_SEGMENT_H

#include <map>
#include <vector>

#include "atom.h"
#include "atomcontainer.h"
#include "qmstate.h"
namespace votca {
namespace xtp {

class Segment : public AtomContainer<Atom> {
 public:
  Segment(std::string name, int id) : AtomContainer<Atom>(name, id){};

  Segment(CheckpointReader& r) : AtomContainer<Atom>("", 0) { ReadFromCpt(r); }

  /// Following notation can be observed in:
  /// [1. Victor, R. et al. Microscopic Simulations of Charge Transport in
  /// Disordered Organic Semiconductors. J. Chem. Theory Comput. 7, 3335–3345
  /// (2011).] Labeling of the following methods follows the following
  /// semantics: U - Energy n - neutral geometry N - neutral state x - excited
  /// geometry X - excited state

  /// UxX - UnN
  void setU_xX_nN(double dU, QMStateType state) {
    _U_xX_nN.setValue(dU, state);
  }
  /// UnX - UnN
  void setU_nX_nN(double dU, QMStateType state) {
    _U_nX_nN.setValue(dU, state);
  }
  /// UxN - UxX
  void setU_xN_xX(double dU, QMStateType state) {
    _U_xN_xX.setValue(dU, state);
  }

  double getU_xX_nN(QMStateType state) const { _U_xX_nN.getValue(state); }

  double getU_nX_nN(QMStateType state) const { _U_nX_nN.getValue(state); }

  double getU_xN_xX(QMStateType state) const { _U_xN_xX.getValue(state); }

  double getSiteEnergy(QMStateType state) const {
    _eMpoles.getValue(state) + _U_xX_nN.getValue(state);
  }

  double getEMpoles(QMStateType state) const { _eMpoles.getValue(state); }

  void setEMpoles(QMStateType state, double energy) {
    _eMpoles.setValue(energy, state);
  }

  void AddMoleculeId(int id) { _molecule_ids.push_back(id); }

  const std::vector<int>& getMoleculeIds() { return _molecule_ids; }

  double getApproxSize() const;

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

 private:
  std::vector<int> _molecule_ids = std::vector<int>(0);

  QMStateCarrierStorage<double> _U_xX_nN;
  QMStateCarrierStorage<double> _U_nX_nN;
  QMStateCarrierStorage<double> _U_xN_xX;
  QMStateCarrierStorage<double> _eMpoles;

  // using caching for approximate size
  mutable double _approxsize = 0.0;
  mutable bool _has_approxsize = false;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENT_H
