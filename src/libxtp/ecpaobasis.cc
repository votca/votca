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

// Standard includes
#include <vector>

// Local VOTCA includes
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/ecpbasisset.h"
#include "votca/xtp/qmmolecule.h"

namespace votca {
namespace xtp {

libecpint::ECP& ECPAOBasis::addShell(const ECPShell& shell,
                                     const QMAtom& atom) {
  libecpint::ECP newU(atom.getPos().data());
  newU.atom_id = int(atom.getId());

  for (const auto& ecpgaussian : shell) {
    newU.addPrimitive(int(ecpgaussian._power), int(shell.getL()), ecpgaussian._decay,
                      ecpgaussian._contraction);
  }
  std::cout<<"shell"<<shell<<std::endl;
  std::cout<<"hi"<<newU<<std::endl;
  _aoshells.push_back(newU);
  std::cout<<"hi2"<<_aoshells.back()<<std::endl;
  return _aoshells.back();
}

Index ECPAOBasis::getMaxL() const {
  int maxL = 0;
  for (const auto& shell : _aoshells) {
    maxL=std::max(shell.L,maxL);
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
  _aoshells.clear();
  _ncore_perAtom.clear();
  _AOBasisSize = 0;
}

void ECPAOBasis::UpdateShellPositions(const QMMolecule& mol) {
  for (libecpint::ECP& shell : _aoshells) {
    const Eigen::Vector3d pos = mol[shell.atom_id].getPos();
    shell.center_ = {pos.x(), pos.y(), pos.z()};
  }
}

void ECPAOBasis::add(const ECPAOBasis& other) {
  Index atomindex_offset = Index(_ncore_perAtom.size());
  for (libecpint::ECP shell : other) {
    shell.atom_id += int(atomindex_offset);
    _AOBasisSize += NumFuncShell(static_cast<L>(shell.L));
    _aoshells.push_back(shell);
  }

  _ncore_perAtom.insert(_ncore_perAtom.end(), other._ncore_perAtom.begin(),
                        other._ncore_perAtom.end());
}

std::vector<std::string> ECPAOBasis::Fill(const ECPBasisSet& bs,
                                          QMMolecule& atoms) {
  _AOBasisSize = 0;
  _aoshells.clear();
  _ncore_perAtom.clear();
  _name = bs.Name();
  std::vector<std::string> non_ecp_elements;
  Index index = 0;
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
      for (const ECPShell& shell : element) {
        addShell(shell, atom);
        index++;
        _AOBasisSize += NumFuncShell(shell.getL());
      }
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

std::ostream& operator<<(std::ostream& out, const libecpint::ECP& shell) {
  out << " Shelltype:" << xtp::EnumToString(static_cast<L>(shell.L))
      << " L:" << Index(shell.L) << " Func:" << shell.N << "\n";
  for (const auto& gaussian : shell.gaussians) {
    out << " Gaussian Decay: " << gaussian.a;
    out << " Power: " << gaussian.n;
    out << " Contractions: " << gaussian.d << "\n";
  }
  return out;
}


std::ostream& operator<<(std::ostream& out, const ECPAOBasis& ecp) {

  out << "Name:" << ecp.Name() << "\n";
  out << " Functions:" << ecp.ECPAOBasisSize()
      << " Shells:" << ecp._aoshells.size() << "\n";
  for (const auto& shell : ecp) {
    out << shell;
  }
  out << "\n"
      << " Atomcharges:";
  for (Index i = 0; i < Index(ecp._ncore_perAtom.size()); i++) {
    out << i << ":" << ecp._ncore_perAtom[i] << " ";
  }
  return out;
}

}  // namespace votca
}  // namespace votca
