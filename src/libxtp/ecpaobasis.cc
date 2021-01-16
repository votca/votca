/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// Standard includes
#include <vector>

// Local VOTCA includes
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/ecpbasisset.h"
#include "votca/xtp/qmmolecule.h"

namespace votca {
namespace xtp {

Index ECPAOBasis::getMaxL() const {
  int maxL = 0;
  for (const auto& potential : _aopotentials) {
    maxL = std::max(potential.getL(), maxL);
  }
  return Index(maxL);
}

void ECPAOBasis::AddECPChargeToMolecule(QMMolecule& mol) const {
  for (Index i = 0; i < mol.size(); i++) {
    mol[i]._ecpcharge = _ncore_perAtom[i];
  }
}

void ECPAOBasis::clear() {
  _name = "";
  _aopotentials.clear();
  _ncore_perAtom.clear();
}

void ECPAOBasis::UpdatePotentialPositions(const QMMolecule& mol) {
  for (libecpint::ECP& potential : _aopotentials) {
    const Eigen::Vector3d pos = mol[potential.atom_id].getPos();
    potential.center_ = {pos.x(), pos.y(), pos.z()};
  }
}

void ECPAOBasis::add(const ECPAOBasis& other) {
  Index atomindex_offset = Index(_ncore_perAtom.size());
  for (libecpint::ECP potential : other) {
    potential.atom_id += int(atomindex_offset);
    _aopotentials.push_back(potential);
  }

  _ncore_perAtom.insert(_ncore_perAtom.end(), other._ncore_perAtom.begin(),
                        other._ncore_perAtom.end());
}

std::vector<std::string> ECPAOBasis::Fill(const ECPBasisSet& bs,
                                          QMMolecule& atoms) {
  _aopotentials.clear();
  _ncore_perAtom.clear();
  _name = bs.Name();
  std::vector<std::string> non_ecp_elements;
  for (QMAtom& atom : atoms) {
    std::string name = atom.getElement();

    bool element_exists = true;

    try {
      bs.getElement(name);
    } catch (std::runtime_error& error) {
      element_exists = false;

      if (std::find(non_ecp_elements.begin(), non_ecp_elements.end(), name) !=
          non_ecp_elements.end()) {
        non_ecp_elements.push_back(name);
      }
    }

    if (element_exists) {
      const ECPElement& element = bs.getElement(name);
      _ncore_perAtom.push_back(element.getNcore());

      libecpint::ECP potential(atom.getPos().data());
      for (const ECPShell& shell : element) {
        for (const auto& gaussian : shell) {
          potential.addPrimitive(int(gaussian._power), int(shell.getL()),
                                 gaussian._decay, gaussian._contraction);
        }
      }
      _aopotentials.push_back(potential);
      _aopotentials.back().atom_id =
          int(atom.getId());  // add atom id here because copyconstructor  of
                              // libecpint::ECP is broken
    } else {
      _ncore_perAtom.push_back(0);
    }
  }

  AddECPChargeToMolecule(atoms);
  return non_ecp_elements;
}

void ECPAOBasis::WriteToCpt(CheckpointWriter& w) const {
  // w(_name, "name");
  // w(_AOBasisSize, "basissize");

  // w(_ncore_perAtom, "atomic ecp charges");
  // Index numofprimitives = 0;
  // for (const auto& shell : _aoshells) {
  //   numofprimitives += shell.getSize();
  // }

  // // this is all to make dummy ECPAOGaussian
  // ECPGaussianPrimitive d(2, 0.1, 0.1);
  // ECPAOGaussianPrimitive dummy2(d);

  // CptTable table = w.openTable("Contractions", dummy2, numofprimitives);

  // std::vector<ECPAOGaussianPrimitive::data> dataVec(numofprimitives);
  // Index i = 0;
  // for (const auto& shell : _aoshells) {
  //   for (const auto& gaussian : shell) {
  //     gaussian.WriteData(dataVec[i], shell);
  //     i++;
  //   }
  // }

  // table.write(dataVec);
}

void ECPAOBasis::ReadFromCpt(CheckpointReader& r) {
  // clear();
  // r(_name, "name");
  // r(_AOBasisSize, "basissize");

  // if (_AOBasisSize > 0) {
  //   r(_ncore_perAtom, "atomic ecp charges");
  //   // this is all to make dummy ECPAOGaussian
  //   ECPGaussianPrimitive d(2, 0.1, 0.1);
  //   ECPAOGaussianPrimitive dummy2(d);

  //   CptTable table = r.openTable("Contractions", dummy2);
  //   std::vector<ECPAOGaussianPrimitive::data> dataVec(table.numRows());
  //   table.read(dataVec);
  //   Index laststartindex = -1;
  //   for (std::size_t i = 0; i < table.numRows(); ++i) {
  //     if (dataVec[i].startindex != laststartindex) {
  //       _aoshells.push_back(ECPAOShell(dataVec[i]));
  //       laststartindex = dataVec[i].startindex;
  //     } else {
  //       _aoshells.back()._gaussians.push_back(
  //           ECPAOGaussianPrimitive(dataVec[i]));
  //     }
  //   }
  // }
}

std::ostream& operator<<(std::ostream& out, const libecpint::ECP& potential) {
  out << " AtomId: " << potential.atom_id << " Components: " << potential.getN()
      << "\n";
  for (const auto& gaussian : potential.gaussians) {
    out << " L: " << gaussian.l;
    out << " Gaussian Decay: " << gaussian.a;
    out << " Power: " << gaussian.n;
    out << " Contractions: " << gaussian.d << "\n";
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const ECPAOBasis& ecp) {

  out << "Name:" << ecp.Name() << "\n";
  out << " Potentials:" << ecp._aopotentials.size() << "\n";
  for (const auto& potential : ecp) {
    out << potential;
  }
  out << "\n"
      << " Atomcharges:";
  for (Index i = 0; i < Index(ecp._ncore_perAtom.size()); i++) {
    out << i << ":" << ecp._ncore_perAtom[i] << " ";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
