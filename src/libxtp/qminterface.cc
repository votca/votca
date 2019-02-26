/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/xtp/qminterface.h>

using boost::format;

namespace votca {
namespace xtp {

ctp::APolarSite *QMInterface::Convert(QMAtom *atm, int id) {

  std::string elem = atm->getType();
  ctp::APolarSite *new_aps = new ctp::APolarSite(id, elem);

  double pol = 0.0;
  try {
    pol = _element.getPolarizability(elem);
  } catch (const std::invalid_argument &out_of_range) {
    std::cout << std::endl
              << "QMMInterface - no default polarizability given "
              << "for element type '" << elem << "'. Defaulting to 1A**3"
              << std::flush;
    pol = 1e-3;
  }
  new_aps->setIsoP(pol);

  tools::vec pos = atm->getPos() * tools::conv::bohr2nm;
  double q = atm->getPartialcharge();
  new_aps->setRank(0);
  new_aps->setPos(pos);
  new_aps->setQ00(q, 0);  // <- charge state 0 <> 'neutral'

  return new_aps;
}

ctp::PolarSeg QMInterface::Convert(std::vector<QMAtom *> &atms) {
  ctp::PolarSeg new_pseg = ctp::PolarSeg();
  std::vector<QMAtom *>::iterator it;
  for (it = atms.begin(); it < atms.end(); ++it) {
    ctp::APolarSite *new_site = this->Convert(*it);
    new_pseg.push_back(new_site);
  }
  return new_pseg;
}

std::vector<QMAtom *> QMInterface::Convert(
    std::vector<ctp::Segment *> segments) {
  std::vector<QMAtom *> qmatoms;
  int AtomId = 0;
  for (ctp::Segment *segment : segments) {
    std::vector<ctp::Atom *> &atoms = segment->Atoms();
    for (ctp::Atom *atom : atoms) {
      if (!atom->HasQMPart()) {
        continue;
      }
      tools::vec pos = atom->getQMPos() * tools::conv::nm2bohr;
      std::string name = atom->getElement();
      QMAtom *qmatom = new QMAtom(AtomId, name, pos);
      AtomId++;
      qmatoms.push_back(qmatom);
    }
  }
  return qmatoms;
}

void QMInterface::GenerateQMAtomsFromPolarSegs(ctp::PolarTop *ptop,
                                               Orbitals &orb) {
  int AtomID = 0;
  // INNER SHELL QM0
  for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
    ctp::PolarSeg *pseg = ptop->QM0()[i];
    for (unsigned int j = 0; j < pseg->size(); ++j) {
      ctp::APolarSite *aps = (*pseg)[j];
      tools::vec pos = aps->getPos() * tools::conv::nm2bohr;
      orb.AddAtom(AtomID, aps->getName(), pos);
      AtomID++;
    }
  }

  return;
}

std::vector<std::shared_ptr<ctp::PolarSeg> > QMInterface::GenerateMultipoleList(
    ctp::PolarTop *ptop) {
  std::vector<std::shared_ptr<ctp::PolarSeg> > MultipoleList;

  // MIDDLE SHELL MM1
  for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
    std::shared_ptr<ctp::PolarSeg> pseg(
        new ctp::PolarSeg(ptop->MM1()[i], false));
    MultipoleList.push_back(pseg);
  }

  // OUTER SHELL MM2
  for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
    std::shared_ptr<ctp::PolarSeg> pseg(
        new ctp::PolarSeg(ptop->MM2()[i], false));
    MultipoleList.push_back(pseg);
  }
  return MultipoleList;
}

void QMInterface::Orbitals2Segment(ctp::Segment &segment,
                                   const Orbitals &orbitals) {

  std::vector<QMAtom *> _atoms = orbitals.QMAtoms();
  std::string type;
  int id = 1;
  for (QMAtom *atom : _atoms) {
    type = atom->getType();
    ctp::Atom *pAtom = new ctp::Atom(id++, type);
    tools::vec position = atom->getPos() * votca::tools::conv::bohr2nm;
    pAtom->setPos(position);
    pAtom->setQMPart(id, position);
    pAtom->setElement(type);
    segment.AddAtom(pAtom);
  }
  return;
}

}  // namespace xtp
}  // namespace votca
