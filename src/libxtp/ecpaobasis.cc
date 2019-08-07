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
#include "votca/xtp/ecpaobasis.h"
#include <vector>
#include <votca/xtp/ecpbasisset.h>
#include <votca/xtp/qmmolecule.h>

namespace votca {
namespace xtp {

ECPAOShell& ECPAOBasis::addShell(const ECPShell& shell, const QMAtom& atom,
                                 int startIndex, int Lmax) {
  _aoshells.push_back(ECPAOShell(shell, atom, startIndex, Lmax));
  return _aoshells.back();
}

std::vector<std::string> ECPAOBasis::Fill(const ECPBasisSet& bs,
                                          QMMolecule& atoms) {
  _AOBasisSize = 0;
  _aoshells.clear();
  _shells_perAtom.clear();
  std::vector<std::vector<int>> shellindex_per_atom;

  std::vector<std::string> non_ecp_elements;
  int index = 0;
  for (QMAtom& atom : atoms) {
    std::string name = atom.getElement();

    bool element_exists = true;
    std::vector<int> shellindex;

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
      atom._ecpcharge = element.getNcore();
      int lmax = element.getLmax();
      for (const ECPShell& shell : element) {
        ECPAOShell& aoshell = addShell(shell, atom, _AOBasisSize, lmax);
        shellindex.push_back(index);
        index++;
        _AOBasisSize += NumFuncShell(shell.getType());
        for (const ECPGaussianPrimitive& gaussian : shell) {
          aoshell.addGaussian(gaussian);
        }
      }
    }
    shellindex_per_atom.push_back(shellindex);
  }

  // have to do it via indeces because if _aoshells resizes pointers are
  // invalidated
  for (const auto& atom_index : shellindex_per_atom) {
    std::vector<const ECPAOShell*> temp;
    for (int index : atom_index) {
      temp.push_back(&_aoshells[index]);
    }
    _shells_perAtom.push_back(temp);
  }

  return non_ecp_elements;
}

}  // namespace xtp
}  // namespace votca
