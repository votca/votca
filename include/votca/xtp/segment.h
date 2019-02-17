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

#include <votca/tools/vec.h>

namespace votca {
namespace xtp {

class Atom;
class Fragment;
class SegmentType;
class Topology;
class Molecule;

class Segment {
 public:
  Segment(int id, std::string name);
  Segment(Segment *stencil);
  ~Segment();

  int getId() const { return _id; }
  const std::string &getName() const { return _name; }

  const tools::vec &getPos() const { return _CoM; }
  void setPos(tools::vec pos) { _CoM = pos; }
  // This gets the center of mass from the MD positions of the atoms
  void calcPos();
  void TranslateBy(const tools::vec &shift);

  void calcApproxSize();
  double getApproxSize() const { return _approxsize; }

  void setHasState(bool yesno, int state);
  bool hasState(int state) const;

  // state: -1 electron +1 hole +2 singlet +3 triplet

  /// Following notation can be observed in:
  /// [1. Victor, R. et al. Microscopic Simulations of Charge Transport in
  /// Disordered Organic Semiconductors. J. Chem. Theory Comput. 7, 3335–3345
  /// (2011).] Labeling of the following methods follows the following
  /// semantics: U - Energy n - neutral geometry N - neutral state x - excited
  /// geometry X - excited state

  /// UxX - UnN
  void setU_xX_nN(double dU, int state);
  /// UnX - UnN
  void setU_nX_nN(double dU, int state);
  /// UxN - UxX
  void setU_xN_xX(double dU, int state);
  double getU_xX_nN(int state) const;
  double getU_nX_nN(int state) const;
  double getU_xN_xX(int state) const;
  double getSiteEnergy(int state) const;

  double getEMpoles(int state) const;
  void setEMpoles(int state, double energy);

  inline void setTopology(Topology *container) { _top = container; }
  Topology *getTopology() { return _top; }
  inline void setMolecule(Molecule *container) { _mol = container; }
  Molecule *getMolecule() { return _mol; }
  inline void setType(SegmentType *type) { _typ = type; }
  SegmentType *getType() { return _typ; }

  void AddFragment(Fragment *fragment);
  void AddAtom(Atom *atom);
  std::vector<Fragment *> &Fragments() { return _fragments; }
  std::vector<Atom *> &Atoms() { return _atoms; }

 private:
  int _id;
  std::string _name;
  SegmentType *_typ;
  Topology *_top;
  Molecule *_mol;

  std::vector<Fragment *> _fragments;
  std::vector<Atom *> _atoms;

  tools::vec _CoM;
  double _approxsize;

  double _U_cC_nN_e;  // from ::EInternal     input     DEFAULT 0
  double _U_cC_nN_h;

  double _U_nC_nN_e;  // from ::EInternal     input     DEFAULT 0
  double _U_nC_nN_h;

  double _U_cN_cC_e;  // from ::EInternal     input     DEFAULT 0
  double _U_cN_cC_h;

  double _U_xX_nN_s;
  double _U_xX_nN_t;

  double _U_nX_nN_s;
  double _U_nX_nN_t;

  double _U_xN_xX_s;
  double _U_xN_xX_t;

  std::vector<double> _eMpoles;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENT_H
