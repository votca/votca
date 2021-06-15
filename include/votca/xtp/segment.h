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
#ifndef VOTCA_XTP_SEGMENT_H
#define VOTCA_XTP_SEGMENT_H

// Standard includes
#include <map>
#include <vector>

// Local VOTCA includes
#include "atom.h"
#include "atomcontainer.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

class Segment : public AtomContainer<Atom> {
 public:
  Segment(std::string name, Index id) : AtomContainer<Atom>(name, id){};
  // cannot use standard AtomContainer constructor because ReadFromCpt is
  // different.
  Segment(CheckpointReader& r) : AtomContainer<Atom>("", 0) { ReadFromCpt(r); }

  ~Segment() override = default;

  /// Following notation can be observed in:
  /// [1. Victor, R. et al. Microscopic Simulations of Charge Transport in
  /// Disordered Organic Semiconductors. J. Chem. Theory Comput. 7, 3335–3345
  /// (2011).] Labeling of the following methods follows the following
  /// semantics: U - Energy n - neutral geometry N - neutral state x - excited
  /// geometry X - excited state

  /// UxX - UnN
  void setU_xX_nN(double dU, QMStateType state) {
    U_xX_nN_.setValue(dU, state);
  }
  /// UnX - UnN
  void setU_nX_nN(double dU, QMStateType state) {
    U_nX_nN_.setValue(dU, state);
  }
  /// UxN - UxX
  void setU_xN_xX(double dU, QMStateType state) {
    U_xN_xX_.setValue(dU, state);
  }

  const Atom* getAtom(Index id) const;

  double getU_xX_nN(QMStateType state) const {
    return U_xX_nN_.getValue(state);
  }

  double getU_nX_nN(QMStateType state) const {
    return U_nX_nN_.getValue(state);
  }

  double getU_xN_xX(QMStateType state) const {
    return U_xN_xX_.getValue(state);
  }

  double getSiteEnergy(QMStateType state) const {
    return site_eng_.getValue(state) + U_xX_nN_.getValue(state);
  }

  double getEMpoles(QMStateType state) const {
    return site_eng_.getValue(state);
  }

  void setEMpoles(QMStateType state, double energy) {
    site_eng_.setValue(energy, state);
  }

  void AddMoleculeId(Index id) { molecule_ids_.push_back(int(id)); }

  const std::vector<Index>& getMoleculeIds() const { return molecule_ids_; }

  double getApproxSize() const;

  void WriteToCpt(CheckpointWriter& w) const override;

  void ReadFromCpt(CheckpointReader& r) override;

  friend std::ostream& operator<<(std::ostream& out, const Segment& container) {
    out << container.getId() << " " << container.getType() << "\n";
    for (const Atom& atom : container) {
      out << atom;
    }
    out << std::endl;
    return out;
  }

 private:
  std::vector<Index> molecule_ids_ = std::vector<Index>(0);

  QMStateCarrierStorage<double> U_xX_nN_;
  QMStateCarrierStorage<double> U_nX_nN_;
  QMStateCarrierStorage<double> U_xN_xX_;
  QMStateCarrierStorage<double> site_eng_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENT_H
