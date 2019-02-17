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

#include <string>
#include <vector>

#include <votca/tools/vec.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/polarsite.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/segmenttype.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace xtp {

/// Default constructor
Segment::Segment(int id, string name)
    : _id(id),
      _name(name),
      _has_e(false),
      _has_h(false),
      _has_s(false),
      _has_t(false) {
  _eMpoles.resize(5);
}

// This constructor creates a copy of the stencil segment, without
// adding it to any containers further up in the hierarchy; i.e. the topology
// and molecules will neither know about the existence of this segment, nor
// be able to access it. Used for creating the ghost in PB corrected pairs.
Segment::Segment(Segment *stencil)
    : _id(stencil->getId()),
      _name(stencil->getName() + "_ghost"),
      _typ(stencil->getType()),
      _top(NULL),
      _mol(NULL),
      _CoM(stencil->getPos()),
      _has_e(false),
      _has_h(false),
      _has_s(false),
      _has_t(false) {
  _eMpoles.resize(5);

  vector<Fragment *>::iterator fit;
  for (fit = stencil->Fragments().begin(); fit < stencil->Fragments().end();
       fit++) {

    Fragment *newFrag = new Fragment(*fit);
    this->AddFragment(newFrag);

    vector<Atom *>::iterator ait;
    for (ait = newFrag->Atoms().begin(); ait < newFrag->Atoms().end(); ait++) {
      this->AddAtom(*ait);
    }
  }
}

Segment::~Segment() {

  vector<Fragment *>::iterator fragit;
  for (fragit = this->Fragments().begin(); fragit < this->Fragments().end();
       fragit++) {
    delete *fragit;
  }
  _fragments.clear();
  _atoms.clear();

  _eMpoles.clear();
}

void Segment::TranslateBy(const vec &shift) {

  _CoM = _CoM + shift;

  vector<Fragment *>::iterator fit;
  for (fit = _fragments.begin(); fit < _fragments.end(); fit++) {

    (*fit)->TranslateBy(shift);
  }
}

void Segment::setHasState(bool yesno, int state) {

  if (state == -1) {
    _has_e = yesno;
  } else if (state == +1) {
    _has_h = yesno;
  } else if (state == +2) {
    _has_s = yesno;
  } else if (state == +3) {
    _has_t = yesno;
  } else {
    throw runtime_error(" ERROR CODE whe__00e11h__");
  }
}

bool Segment::hasState(int state) const {
  bool result;
  if (state == -1) {
    result = _has_e;
  } else if (state == +1) {
    result = _has_h;
  } else if (state == +2) {
    result = _has_s;
  } else if (state == +3) {
    result = _has_t;
  } else {
    throw runtime_error(" ERROR CODE whe__00s11o__");
  }
  return result;
}

void Segment::setOcc(double occ, int e_h_s_t) {

  if (e_h_s_t == -1) {
    _occ_e = occ;
  } else if (e_h_s_t == +1) {
    _occ_h = occ;
  } else if (e_h_s_t == +2) {
    _occ_s = occ;
  } else if (e_h_s_t == +3) {
    _occ_t = occ;
  } else {
    throw runtime_error(" ERROR CODE whe__00s11o__");
  }
}

double Segment::getOcc(int e_h_s_t) const {
  double result;
  if (e_h_s_t == -1) {
    result = _occ_e;
  } else if (e_h_s_t == +1) {
    result = _occ_h;
  } else if (e_h_s_t == +2) {
    result = _occ_s;
  } else if (e_h_s_t == +3) {
    result = _occ_t;
  } else {
    throw runtime_error(" ERROR CODE whe__00s11o__");  // blabla what do I do
                                                       // here?
  }
  return result;
}

void Segment::setU_cC_nN(double dU, int state) {

  if (state == -1) {
    _U_cC_nN_e = dU;
  } else if (state == +1) {
    _U_cC_nN_h = dU;
  } else {
    throw runtime_error(" ERROR CODE whe__00u11a__");
  }
}

void Segment::setU_nC_nN(double dU, int state) {

  if (state == -1) {
    _U_nC_nN_e = dU;
  } else if (state == +1) {
    _U_nC_nN_h = dU;
  } else {
    throw runtime_error(" ERROR CODE whe__00u11b__");
  }
}

void Segment::setU_cN_cC(double dU, int state) {

  if (state == -1) {
    _U_cN_cC_e = dU;
  } else if (state == +1) {
    _U_cN_cC_h = dU;
  } else {
    throw runtime_error(" ERROR CODE whe__00u11c__");
  }
}

void Segment::setU_xX_nN(double dU, int state) {

  if (state == +2) {
    _U_xX_nN_s = dU;
  } else if (state == +3) {
    _U_xX_nN_t = dU;
  } else {
    throw runtime_error(" ERROR CODE whe__00u11d__");  // blabla?? What do I do
                                                       // here?
  }
}

void Segment::setU_nX_nN(double dU, int state) {

  if (state == +2) {
    _U_nX_nN_s = dU;
  } else if (state == +3) {
    _U_nX_nN_t = dU;
  } else {
    throw runtime_error(" ERROR CODE whe__00u11d__");  // blabla?? What do I do
                                                       // here?
  }
}

void Segment::setU_xN_xX(double dU, int state) {

  if (state == +2) {
    _U_xN_xX_s = dU;
  } else if (state == +3) {
    _U_xN_xX_t = dU;
  } else {
    throw runtime_error(" ERROR CODE whe__00u11d__");  // blabla?? What do I do
                                                       // here?
  }
}

double Segment::getU_xX_nN(int state) const {

  return (state == +3) ? _U_xX_nN_t : _U_xX_nN_s;
}

double Segment::getU_nX_nN(int state) const {

  return (state == +3) ? _U_nX_nN_t : _U_nX_nN_s;
}

double Segment::getU_xN_xX(int state) const {

  return (state == +3) ? _U_xN_xX_t : _U_xN_xX_s;
}

double Segment::getU_cC_nN(int state) const {

  return (state == -1) ? _U_cC_nN_e : _U_cC_nN_h;
}

double Segment::getU_nC_nN(int state) const {

  return (state == -1) ? _U_nC_nN_e : _U_nC_nN_h;
}

double Segment::getU_cN_cC(int state) const {

  return (state == -1) ? _U_cN_cC_e : _U_cN_cC_h;
}

double Segment::getSiteEnergy(int state) const {

  double result;
  if (state == -1) {
    result = getEMpoles(state) + _U_cC_nN_e;
  } else if (state == +1) {
    result = getEMpoles(state) + _U_cC_nN_h;
  } else if (state == +2) {
    result = getEMpoles(state) + _U_xX_nN_s;
  } else if (state == +3) {
    result = getEMpoles(state) + _U_xX_nN_t;
  } else {
    throw runtime_error(" ERROR CODE whe__00s11o__");  // blabla what do I do
                                                       // here?
  }
  return result;
}

void Segment::setEMpoles(int state, double energy) {

  _hasChrgState.resize(5);
  _hasChrgState[state + 1] = true;
  _eMpoles[state + 1] = energy;
}

double Segment::getEMpoles(int state) const {
  return _eMpoles[state + 1] - _eMpoles[1];
}

void Segment::AddFragment(Fragment *fragment) {
  _fragments.push_back(fragment);
  fragment->setSegment(this);
}

void Segment::AddAtom(Atom *atom) {
  _atoms.push_back(atom);
  atom->setSegment(this);
}

void Segment::calcPos() {
  vec pos = vec(0, 0, 0);
  double totWeight = 0.0;

  for (unsigned int i = 0; i < _atoms.size(); i++) {
    pos += _atoms[i]->getPos() * _atoms[i]->getWeight();
    totWeight += _atoms[i]->getWeight();
  }

  _CoM = pos / totWeight;
}

void Segment::calcApproxSize() {
  _approxsize = 0.0;

  tools::vec min = vec(numeric_limits<double>::max());
  tools::vec max = vec(numeric_limits<double>::min());
  vector<Fragment *>::iterator fragit1;
  for (fragit1 = Fragments().begin(); fragit1 < Fragments().end(); fragit1++) {
    const tools::vec &pos = (*fragit1)->getPos();
    if (pos.getX() > max.getX()) {
      max.x() = pos.getX();
    } else if (pos.getX() < min.getX()) {
      min.x() = pos.getX();
    }
    if (pos.getY() > max.getY()) {
      max.y() = pos.getY();
    } else if (pos.getY() < min.getY()) {
      min.y() = pos.getY();
    }
    if (pos.getZ() > max.getZ()) {
      max.z() = pos.getZ();
    } else if (pos.getZ() < min.getZ()) {
      min.z() = pos.getZ();
    }
  }

  _approxsize = abs(max - min);

  return;
}

}  // namespace xtp
}  // namespace votca
