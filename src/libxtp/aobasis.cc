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
#include "votca/xtp/aobasis.h"
#include <vector>
#include <votca/xtp/basisset.h>
#include <votca/xtp/qmmolecule.h>

namespace votca {
namespace xtp {

AOShell& AOBasis::addShell(const Shell& shell, const QMAtom& atom,
                           Index startIndex) {
  _aoshells.push_back(AOShell(shell, atom, startIndex));
  return _aoshells.back();
}

void AOBasis::ReorderMOs(Eigen::MatrixXd& v, const std::string& start,
                         const std::string& target) const {

  if (start == target) {
    return;
  }

  if (target == "orca" || target == "nwchem") {
    std::vector<Index> multiplier = getMultiplierVector(target, start);
    // and reorder rows of _orbitals->_mo_coefficients() accordingly
    MultiplyMOs(v, multiplier);
  }

  // get reordering vector _start -> target
  std::vector<Index> order = getReorderVector(start, target);

  // Sanity check
  if (v.rows() != Index(order.size())) {
    throw std::runtime_error("Size mismatch in ReorderMOs " +
                             std::to_string(v.rows()) + ":" +
                             std::to_string(order.size()));
  }

  // actual swapping of coefficients
  for (Index s = 1, d; s < (Index)order.size(); ++s) {
    for (d = order[s]; d < s; d = order[d]) {
      ;
    }
    if (d == s) {
      while (d = order[d], d != s) {
        v.row(s).swap(v.row(d));
      }
    }
  }
  // NWChem has some strange minus in d-functions
  if (start == "nwchem" || start == "orca") {
    std::vector<Index> multiplier = getMultiplierVector(start, target);
    // and reorder rows of _orbitals->_mo_coefficients() accordingly
    MultiplyMOs(v, multiplier);
  }
  return;
}

void AOBasis::MultiplyMOs(Eigen::MatrixXd& v,
                          const std::vector<Index>& multiplier) const {
  // Sanity check
  if (v.cols() != Index(multiplier.size())) {
    std::cerr << "Size mismatch in MultiplyMOs" << v.cols() << ":"
              << multiplier.size() << std::endl;
    throw std::runtime_error("Abort!");
  }
  for (Index i_basis = 0; i_basis < v.cols(); i_basis++) {
    v.row(i_basis) = multiplier[i_basis] * v.row(i_basis);
  }
  return;
}

std::vector<Index> AOBasis::getMultiplierVector(
    const std::string& start, const std::string& target) const {
  std::vector<Index> multiplier;
  multiplier.reserve(_AOBasisSize);
  std::string s = start;
  std::string t = target;
  if (start == "xtp") {
    s = target;
    t = start;
  }
  // go through basisset
  for (const AOShell& shell : (*this)) {
    addMultiplierShell(s, t, shell.getType(), multiplier);
  }
  return multiplier;
}

void AOBasis::addMultiplierShell(const std::string& start,
                                 const std::string& target,
                                 const std::string& shell_type,
                                 std::vector<Index>& multiplier) const {
  // multipliers were all found using code, hard to establish

  if (target == "xtp") {
    // single type shells defined here
    if (shell_type.length() == 1) {
      if (shell_type == "S") {
        multiplier.push_back(1);
      } else if (shell_type == "P") {
        multiplier.push_back(1);
        multiplier.push_back(1);
        multiplier.push_back(1);
      } else if (shell_type == "D") {
        if (start == "nwchem") {
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(-1);
          multiplier.push_back(1);
          multiplier.push_back(1);
        } else if (start == "orca") {
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
        } else {
          std::cerr << "Tried to get multipliers d-functions from package "
                    << start << ".";
          throw std::runtime_error("Multiplication not implemented yet!");
        }
      } else if (shell_type == "F") {
        if (start == "orca") {
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(-1);
          multiplier.push_back(-1);

        } else if (start == "nwchem") {
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(-1);
          multiplier.push_back(+1);
          multiplier.push_back(+1);
          multiplier.push_back(+1);
          multiplier.push_back(-1);
        } else {
          std::cerr << "Tried to get multipliers f-functions from package "
                    << start << ".";
          throw std::runtime_error("Multiplication not implemented yet!");
        }
      } else if (shell_type == "G") {
        if (start == "orca") {
          // Not checked yet
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(1);
          multiplier.push_back(-1);
          multiplier.push_back(-1);
          multiplier.push_back(-1);
          multiplier.push_back(-1);

        } else {
          std::cerr << "Tried to get multipliers g-functions from package "
                    << start << ".";
          throw std::runtime_error("Multiplication not implemented yet!");
        }
      } else {
        std::cerr << "Tried to get multipliers h-functions . ";
        throw std::runtime_error("Multiplication not implemented yet!");
      }
    } else {
      // for combined shells, iterate over all contributions
      for (Index i = 0; i < Index(shell_type.length()); ++i) {
        std::string local_shell = std::string(shell_type, i, 1);
        addMultiplierShell(start, target, local_shell, multiplier);
      }
    }
  } else {

    std::cerr << "Tried to reorder functions (multiplier)  from " << start
              << " to " << target << std::endl;
    throw std::runtime_error("Reordering not implemented yet!");
  }
  return;
}

std::vector<Index> AOBasis::getReorderVector(const std::string& start,
                                             const std::string& target) const {
  std::vector<Index> neworder;
  neworder.reserve(_AOBasisSize);
  std::string s;
  std::string t;
  if (start == "xtp") {
    s = target;
    t = start;
  } else {
    s = start;
    t = target;
  }
  // go through basisset
  for (const AOShell& shell : (*this)) {
    addReorderShell(s, t, shell.getType(), neworder);
  }
  if (start == "xtp") {
    neworder = invertOrder(neworder);
  }
  return neworder;
}

std::vector<Index> AOBasis::invertOrder(const std::vector<Index>& order) const {

  std::vector<Index> neworder = std::vector<Index>(order.size());
  for (Index i = 0; i < Index(order.size()); i++) {
    neworder[order[i]] = Index(i);
  }
  return neworder;
}

void AOBasis::addReorderShell(const std::string& start,
                              const std::string& target,
                              const std::string& shell_type,
                              std::vector<Index>& order) const {
  // Reordering is given by email from gaussian, orca output MOs, and
  // http://www.nwchem-sw.org/index.php/Release66:Basis for nwchem

  // current length of vector

  Index cur_pos = Index(order.size()) - 1;

  if (target == "xtp") {
    // single type shells defined here
    if (shell_type.length() == 1) {
      if (shell_type == "S") {
        order.push_back(cur_pos + 1);
      }  // for S

      // votca order is z,y,x e.g. Y1,0 Y1,-1 Y1,1
      else if (shell_type == "P") {
        if (start == "orca") {
          // orca order is z,x,y Y1,0,Y1,1,Y1,-1
          order.push_back(cur_pos + 1);
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 2);
        } else if (start == "gaussian" || start == "nwchem") {
          // nwchem gaussian x,y,z Y1,1 Y1,-1 Y1,0
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 2);
          order.push_back(cur_pos + 1);
        } else {
          std::cerr << "Tried to reorder p-functions from package " << start
                    << ".";
          throw std::runtime_error("Reordering not implemented yet!");
        }
      }  // for P
         // votca order is d3z2-r2 dyz dxz dxy dx2-y2 e.g. Y2,0 Y2,-1 Y2,1 Y2,-2
         // Y2,2
      else if (shell_type == "D") {
        // orca order is d3z2-r2 dxz dyz dx2-y2 dxy e.g. Y2,0 Y2,1 Y2,-1 Y2,2
        // Y2,-2
        if (start == "gaussian" || start == "orca") {
          order.push_back(cur_pos + 1);
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 2);
          order.push_back(cur_pos + 5);
          order.push_back(cur_pos + 4);
        } else if (start == "nwchem") {
          // nwchem order is dxy dyz d3z2-r2 -dxz dx2-y2, e.g. Y2,-2 Y2,-1 Y2,0
          // Y2,1 Y2,2
          order.push_back(cur_pos + 4);
          order.push_back(cur_pos + 2);
          order.push_back(cur_pos + 1);
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 5);
        } else {
          std::cerr << "Tried to reorder d-functions from package " << start
                    << ".";
          throw std::runtime_error("Reordering not implemented yet!");
        }
      } else if (shell_type == "F") {
        // ordering for votca is Yl,0 Yl,-1 Yl,1 ......Yl,-m Yl,m
        if (start == "gaussian" || start == "orca") {
          // ordering for gaussian and orca is Yl,0 Yl,1 Yl,-1 ......Yl,m Yl,-m
          order.push_back(cur_pos + 1);
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 2);
          order.push_back(cur_pos + 5);
          order.push_back(cur_pos + 4);
          order.push_back(cur_pos + 7);
          order.push_back(cur_pos + 6);
        } else if (start == "nwchem") {
          // ordering for nwchem is fxxy-yyy,
          // fxyz,fyzz-xxy-yyy,fzzz-xxz-yyz,f-xzz+xxx+xyy,fxxz-yyz,fxyy-xxx
          // e.g. Y3,-3 Y3,-2 Y3,-1 Y3,0 Y3,1 Y3,2 Y3,3
          order.push_back(cur_pos + 6);
          order.push_back(cur_pos + 4);
          order.push_back(cur_pos + 2);
          order.push_back(cur_pos + 1);
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 5);
          order.push_back(cur_pos + 7);
        } else {
          std::cerr << "Tried to reorder f-functions from package " << start
                    << ".";
          throw std::runtime_error("Reordering not implemented yet!");
        }
      } else if (shell_type == "G") {
        // ordering for votca is Yl,0 Yl,-1 Yl,1 ......Yl,-m Yl,m
        if (start == "gaussian" || start == "orca") {
          // ordering for gaussian and orca is Yl,0 Yl,1 Yl,-1 ......Yl,m Yl,-m
          order.push_back(cur_pos + 1);
          order.push_back(cur_pos + 3);
          order.push_back(cur_pos + 2);
          order.push_back(cur_pos + 5);
          order.push_back(cur_pos + 4);
          order.push_back(cur_pos + 7);
          order.push_back(cur_pos + 6);
          order.push_back(cur_pos + 9);
          order.push_back(cur_pos + 8);
        } else {
          std::cerr << "Tried to reorder g-functions from package " << start
                    << ".";
          throw std::runtime_error("Reordering not implemented");
        }
      } else {
        std::cerr << "Tried to reorder functions  of shell type " << shell_type
                  << std::endl;
        throw std::runtime_error("Reordering not implemented");
      }
    } else {
      // for combined shells, iterate over all contributions
      //_nbf = 0;
      for (Index i = 0; i < Index(shell_type.length()); ++i) {
        std::string local_shell = std::string(shell_type, i, 1);
        this->addReorderShell(start, target, local_shell, order);
      }
    }

  } else {
    std::cerr << "Tried to reorder functions (neworder) from " << start
              << " to " << target << std::endl;
    throw std::runtime_error("Reordering not implemented yet!");
  }
  return;
}

const std::vector<const AOShell*> AOBasis::getShellsofAtom(Index AtomId) const {
  std::vector<const AOShell*> result;
  for (const auto& aoshell : _aoshells) {
    if (aoshell.getAtomIndex() == AtomId) {
      result.push_back(&aoshell);
    }
  }
  return result;
}

void AOBasis::Fill(const BasisSet& bs, const QMMolecule& atoms) {
  _AOBasisSize = 0;
  _aoshells.clear();
  _FuncperAtom.clear();
  // loop over atoms
  for (const QMAtom& atom : atoms) {
    Index atomfunc = 0;
    const std::string& name = atom.getElement();
    const Element& element = bs.getElement(name);
    for (const Shell& shell : element) {
      Index numfuncshell = NumFuncShell(shell.getType());
      AOShell& aoshell = addShell(shell, atom, _AOBasisSize);
      _AOBasisSize += numfuncshell;
      atomfunc += numfuncshell;
      for (const GaussianPrimitive& gaussian : shell) {
        aoshell.addGaussian(gaussian);
      }
      aoshell.CalcMinDecay();
      aoshell.normalizeContraction();
    }
    _FuncperAtom.push_back(atomfunc);
  }
  return;
}

}  // namespace xtp
}  // namespace votca
