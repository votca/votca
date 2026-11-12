/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#include <algorithm>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/activeregion.h"
#include "votca/xtp/aomatrix.h"

namespace votca {
namespace xtp {

std::vector<std::vector<Index>> ActiveRegion::BondGraph(const QMMolecule& mol) {
  tools::Elements elements;
  const Index n = mol.size();
  std::vector<double> radius(static_cast<std::size_t>(n));
  for (Index i = 0; i < n; ++i) {
    radius[std::size_t(i)] = elements.getCovRad(mol[i].getElement(), "bohr");
  }
  std::vector<std::vector<Index>> graph(static_cast<std::size_t>(n));
  for (Index i = 0; i < n; ++i) {
    for (Index j = i + 1; j < n; ++j) {
      const double d = (mol[i].getPos() - mol[j].getPos()).norm();
      if (d < 1.2 * (radius[std::size_t(i)] + radius[std::size_t(j)])) {
        graph[std::size_t(i)].push_back(j);
        graph[std::size_t(j)].push_back(i);
      }
    }
  }
  return graph;
}

Eigen::MatrixXd ActiveRegion::AtomPopulations(const Orbitals& orb,
                                              const Eigen::MatrixXd& mos) {
  const AOBasis& basis = orb.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  // Mulliken: q_mu = c_mu (S c)_mu, summed over the functions of an atom
  const Eigen::MatrixXd q = mos.cwiseProduct(overlap.Matrix() * mos);
  const std::vector<Index>& per_atom = basis.getFuncPerAtom();
  Eigen::MatrixXd pop(Index(per_atom.size()), mos.cols());
  Index start = 0;
  for (std::size_t a = 0; a < per_atom.size(); ++a) {
    pop.row(Index(a)) = q.middleRows(start, per_atom[a]).colwise().sum();
    start += per_atom[a];
  }
  return pop;
}

namespace {
// connected components of the subgraph on `nodes`
std::vector<std::vector<Index>> Components(
    const std::vector<std::vector<Index>>& graph,
    const std::vector<bool>& in_set) {
  const Index n = Index(graph.size());
  std::vector<Index> label(std::size_t(n), -1);
  std::vector<std::vector<Index>> comps;
  for (Index s = 0; s < n; ++s) {
    if (!in_set[std::size_t(s)] || label[std::size_t(s)] >= 0) {
      continue;
    }
    comps.emplace_back();
    std::vector<Index> stack{s};
    label[std::size_t(s)] = Index(comps.size()) - 1;
    while (!stack.empty()) {
      const Index i = stack.back();
      stack.pop_back();
      comps.back().push_back(i);
      for (Index j : graph[std::size_t(i)]) {
        if (in_set[std::size_t(j)] && label[std::size_t(j)] < 0) {
          label[std::size_t(j)] = label[std::size_t(s)];
          stack.push_back(j);
        }
      }
    }
  }
  return comps;
}
}  // namespace

void ActiveRegion::CheckConnected(const QMMolecule& mol,
                                  const std::vector<Index>& active) {
  const auto graph = BondGraph(mol);
  const Index n = mol.size();
  std::vector<bool> all(std::size_t(n), true);
  std::vector<bool> act(std::size_t(n), false);
  for (Index a : active) {
    if (a < 0 || a >= n) {
      throw std::runtime_error("dft_in_dft: active atom index " +
                               std::to_string(a) + " out of range");
    }
    act[std::size_t(a)] = true;
  }
  for (const auto& molecule : Components(graph, all)) {
    std::vector<bool> act_here(std::size_t(n), false);
    for (Index i : molecule) {
      act_here[std::size_t(i)] = act[std::size_t(i)];
    }
    const auto pieces = Components(graph, act_here);
    if (pieces.size() > 1) {
      std::ostringstream msg;
      msg << "dft_in_dft: the active atoms split a molecule into "
          << pieces.size() << " disconnected pieces:";
      for (const auto& p : pieces) {
        msg << " {";
        for (std::size_t k = 0; k < p.size(); ++k) {
          msg << (k ? "," : "") << p[k];
        }
        msg << "}";
      }
      msg << ". Each molecule should contribute one connected piece.";
      throw std::runtime_error(msg.str());
    }
  }
}

ActiveRegion::Report ActiveRegion::Evaluate(const Orbitals& full,
                                            const std::vector<Index>& active,
                                            const Options& options) {
  const QMMolecule& mol = full.QMAtoms();
  const Index n = mol.size();
  CheckConnected(mol, active);

  // frontier orbitals of the full calculation
  const Index homo = full.getNumberOfAlphaElectrons() - 1;
  const Eigen::VectorXd& e = full.MOs().eigenvalues();
  std::vector<Index> orbs;
  for (Index k = homo; k >= 0 && e(k) >= e(homo) - options.window; --k) {
    orbs.push_back(k);
  }
  if (options.include_virtuals) {
    const Index lumo = homo + 1;
    for (Index k = lumo; k < e.size() && e(k) <= e(lumo) + options.window;
         ++k) {
      orbs.push_back(k);
    }
  }
  Eigen::MatrixXd mos(full.MOs().eigenvectors().rows(), Index(orbs.size()));
  Report report;
  report.energies.resize(Index(orbs.size()));
  for (std::size_t k = 0; k < orbs.size(); ++k) {
    mos.col(Index(k)) = full.MOs().eigenvectors().col(orbs[k]);
    report.energies(Index(k)) = e(orbs[k]);
  }
  report.orbitals = orbs;
  const Eigen::MatrixXd pop = AtomPopulations(full, mos);
  // the largest population of an atom over the frontier orbitals
  const Eigen::VectorXd maxpop = pop.rowwise().maxCoeff();

  std::vector<bool> act(std::size_t(n), false);
  for (Index a : active) {
    act[std::size_t(a)] = true;
  }
  auto coverage = [&](const std::vector<bool>& in) {
    Eigen::VectorXd c = Eigen::VectorXd::Zero(pop.cols());
    for (Index i = 0; i < n; ++i) {
      if (in[std::size_t(i)]) {
        c += pop.row(i).transpose();
      }
    }
    return c;
  };
  report.coverage_input = coverage(act);

  const auto graph = BondGraph(mol);
  if (options.grow) {
    const std::vector<bool> input = act;
    // every atom above the threshold in a molecule that has active atoms
    // joins, together with the shortest bond path that connects it to the
    // region (an atom on a node of the orbital can lie on that path)
    std::vector<bool> molecule_active(std::size_t(n), false);
    {
      std::vector<bool> all(std::size_t(n), true);
      for (const auto& molecule : Components(graph, all)) {
        bool any = false;
        for (Index i : molecule) {
          any = any || act[std::size_t(i)];
        }
        for (Index i : molecule) {
          molecule_active[std::size_t(i)] = any;
        }
      }
    }
    bool changed = true;
    while (changed) {
      changed = false;
      // breadth-first search from the region to the nearest atom above
      // the threshold that is not yet active
      std::vector<Index> parent(std::size_t(n), -2);
      std::vector<Index> queue;
      for (Index i = 0; i < n; ++i) {
        if (act[std::size_t(i)]) {
          parent[std::size_t(i)] = -1;
          queue.push_back(i);
        }
      }
      Index target = -1;
      for (std::size_t q = 0; q < queue.size() && target < 0; ++q) {
        for (Index j : graph[std::size_t(queue[q])]) {
          if (parent[std::size_t(j)] != -2) {
            continue;
          }
          parent[std::size_t(j)] = queue[q];
          if (maxpop(j) > options.threshold &&
              molecule_active[std::size_t(j)]) {
            target = j;
            break;
          }
          queue.push_back(j);
        }
      }
      for (Index i = target; i >= 0 && !act[std::size_t(i)];
           i = parent[std::size_t(i)]) {
        act[std::size_t(i)] = true;
        changed = true;
      }
    }
    // heavy atoms bonded to two or more active atoms (a ring atom between
    // two active neighbours, an atom enclosed by the region)
    changed = true;
    while (changed) {
      changed = false;
      for (Index i = 0; i < n; ++i) {
        if (act[std::size_t(i)] || mol[i].getElement() == "H") {
          continue;
        }
        Index active_neighbours = 0;
        for (Index j : graph[std::size_t(i)]) {
          active_neighbours += act[std::size_t(j)] ? 1 : 0;
        }
        if (active_neighbours >= 2) {
          act[std::size_t(i)] = true;
          changed = true;
        }
      }
    }
    // hydrogens of active heavy atoms (also of the input ones: an X-H bond
    // should never be cut)
    for (Index i = 0; i < n; ++i) {
      if (act[std::size_t(i)] && mol[i].getElement() != "H") {
        for (Index j : graph[std::size_t(i)]) {
          if (mol[j].getElement() == "H") {
            act[std::size_t(j)] = true;
          }
        }
      }
    }
    for (Index i = 0; i < n; ++i) {
      if (act[std::size_t(i)] && !input[std::size_t(i)]) {
        report.added.push_back(i);
      }
    }
  }
  for (Index i = 0; i < n; ++i) {
    if (act[std::size_t(i)]) {
      report.atoms.push_back(i);
    } else if (maxpop(i) > options.threshold) {
      report.unreachable.push_back(i);
    }
  }
  if (!options.grow) {
    report.unreachable.clear();  // only meaningful after growing
  }
  for (Index i = 0; i < n; ++i) {
    if (!act[std::size_t(i)] && mol[i].getElement() == "H") {
      for (Index j : graph[std::size_t(i)]) {
        if (act[std::size_t(j)]) {
          report.cut_hydrogens.push_back(i);
          break;
        }
      }
    }
  }
  for (Index i = 0; i < n; ++i) {
    if (act[std::size_t(i)] || mol[i].getElement() == "H") {
      continue;
    }
    Index active_neighbours = 0;
    for (Index j : graph[std::size_t(i)]) {
      active_neighbours += act[std::size_t(j)] ? 1 : 0;
    }
    if (active_neighbours >= 2) {
      report.enclosed.push_back(i);
    }
  }
  report.coverage_final = coverage(act);
  CheckConnected(mol, report.atoms);
  return report;
}

std::string ActiveRegion::Report::Format() const {
  std::ostringstream out;
  out << "Frontier orbitals of the full calculation on the active region "
         "(Mulliken population):\n";
  for (std::size_t k = 0; k < orbitals.size(); ++k) {
    out << "  MO " << std::setw(5) << orbitals[k] << "  " << std::fixed
        << std::setprecision(4) << energies(Index(k)) * tools::conv::hrt2ev
        << " eV  input " << std::setprecision(3) << coverage_input(Index(k));
    if (!added.empty()) {
      out << "  grown " << coverage_final(Index(k));
    }
    out << "\n";
  }
  if (!added.empty()) {
    out << "  added " << added.size() << " atoms:";
    for (Index a : added) {
      out << " " << a;
    }
    out << "\n";
  }
  if (!cut_hydrogens.empty()) {
    out << "  warning: inactive hydrogens bonded to active atoms (X-H bonds "
           "cut; grow adds them):";
    for (Index a : cut_hydrogens) {
      out << " " << a;
    }
    out << "\n";
  }
  if (!enclosed.empty()) {
    out << "  warning: inactive atoms bonded to two or more active atoms "
           "(several bonds cut):";
    for (Index a : enclosed) {
      out << " " << a;
    }
    out << "\n";
  }
  if (!unreachable.empty()) {
    out << "  above the threshold in molecules without active atoms:";
    for (Index a : unreachable) {
      out << " " << a;
    }
    out << "\n";
  }
  return out.str();
}

}  // namespace xtp
}  // namespace votca
