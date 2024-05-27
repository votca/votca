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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 *  olee_ the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Third party includes
#include <boost/math/constants/constants.hpp>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/radial_euler_maclaurin_rule.h"

namespace votca {
namespace xtp {

std::vector<double> EulerMaclaurinGrid::CalculatePruningIntervals(
    const std::string& element) {
  std::vector<double> r;
  // get Bragg-Slater Radius for this element
  double BSradius = BraggSlaterRadii_.at(element);
  // row type of element
  Index RowType = pruning_set_.at(element);

  if (RowType == 1) {
    r.push_back(0.25 * BSradius);
    r.push_back(0.5 * BSradius);
    r.push_back(1.0 * BSradius);
    r.push_back(4.5 * BSradius);
  } else if (RowType == 2) {
    r.push_back(0.1667 * BSradius);
    r.push_back(0.5 * BSradius);
    r.push_back(0.9 * BSradius);
    r.push_back(3.5 * BSradius);
  } else if (RowType == 3) {
    r.push_back(0.1 * BSradius);
    r.push_back(0.4 * BSradius);
    r.push_back(0.8 * BSradius);
    r.push_back(2.5 * BSradius);
  } else {
    throw std::runtime_error(
        "EulerMaclaurinGrid::CalculatePruningIntervals:Pruning unsupported for "
        "RowType");
  }
  return r;
}

void EulerMaclaurinGrid::FillElementRangeMap(const AOBasis& aobasis,
                                             const QMMolecule& atoms,
                                             double eps) {
  std::map<std::string, min_exp>::iterator it;
  for (const QMAtom& atom : atoms) {
    std::string name = atom.getElement();
    // is this element already in map?
    it = element_ranges_.find(name);
    // only proceed, if element data does not exist yet
    if (it == element_ranges_.end()) {
      min_exp this_atom;
      double range_max = std::numeric_limits<double>::min();
      double decaymin = std::numeric_limits<double>::max();
      Index lvalue = std::numeric_limits<Index>::min();
      const std::vector<const AOShell*> shells =
          aobasis.getShellsofAtom(atom.getId());
      // and loop over all shells to figure out minimum decay constant and
      // angular momentum of this function
      for (const AOShell* shell : shells) {
        Index lmax = Index(shell->getL());
        if (shell->getMinDecay() < decaymin) {
          decaymin = shell->getMinDecay();
          lvalue = lmax;
        }
        double range = DetermineCutoff(2 * decaymin, 2 * lvalue + 2, eps);
        if (range > range_max) {
          this_atom.alpha = decaymin;
          this_atom.l = lvalue;
          this_atom.range = range;
          range_max = range;
        }
      }  // shells
      element_ranges_[name] = this_atom;
    }  // new element
  }  // atoms
}

void EulerMaclaurinGrid::RefineElementRangeMap(const AOBasis& aobasis,
                                               const QMMolecule& atoms,
                                               double eps) {
  AOOverlap overlap;
  overlap.Fill(aobasis);

  // get collapsed index list
  std::vector<Index> idxstart;
  const std::vector<Index>& idxsize = aobasis.getFuncPerAtom();
  Index start = 0;
  for (Index size : idxsize) {
    idxstart.push_back(start);
    start += size;
  }
  // refining by going through all atom combinations
  for (Index i = 0; i < atoms.size(); ++i) {
    const QMAtom& atom_a = atoms[i];
    Index a_start = idxstart[i];
    Index a_size = idxsize[i];
    double range_max = std::numeric_limits<double>::min();
    // get preset values for this atom type
    double alpha_a = element_ranges_.at(atom_a.getElement()).alpha;
    Index l_a = element_ranges_.at(atom_a.getElement()).l;
    const Eigen::Vector3d& pos_a = atom_a.getPos();
    // Cannot iterate only over j<i because it is not symmetric due to shift_2g
    for (Index j = 0; j < atoms.size(); ++j) {
      if (i == j) {
        continue;
      }
      const QMAtom& atom_b = atoms[j];
      Index b_start = idxstart[j];
      Index b_size = idxsize[j];
      const Eigen::Vector3d& pos_b = atom_b.getPos();
      // find overlap block of these two atoms
      Eigen::MatrixXd overlapblock =
          overlap.Matrix().block(a_start, b_start, a_size, b_size);
      // determine abs max of this block
      double s_max = overlapblock.cwiseAbs().maxCoeff();

      if (s_max > 1e-5) {
        double range = DetermineCutoff(
            alpha_a + element_ranges_.at(atom_b.getElement()).alpha,
            l_a + element_ranges_.at(atom_b.getElement()).l + 2, eps);
        // now do some update trickery from Gaussian product formula
        double dist = (pos_b - pos_a).norm();
        double shift_2g =
            dist * alpha_a /
            (alpha_a + element_ranges_.at(atom_b.getElement()).alpha);
        range += (shift_2g + dist);
        if (range > range_max) {
          range_max = range;
        }
      }
    }
    if (std::round(range_max) > element_ranges_.at(atom_a.getElement()).range) {
      element_ranges_.at(atom_a.getElement()).range = std::round(range_max);
    }
  }
}

void EulerMaclaurinGrid::CalculateRadialCutoffs(const AOBasis& aobasis,
                                                const QMMolecule& atoms,
                                                const std::string& gridtype) {

  double eps = Accuracy[gridtype];
  FillElementRangeMap(aobasis, atoms, eps);
  RefineElementRangeMap(aobasis, atoms, eps);
  return;
}

std::map<std::string, GridContainers::radial_grid>
    EulerMaclaurinGrid::CalculateAtomicRadialGrids(const AOBasis& aobasis,
                                                   const QMMolecule& atoms,
                                                   const std::string& type) {

  CalculateRadialCutoffs(aobasis, atoms, type);
  std::map<std::string, GridContainers::radial_grid> result;
  for (const auto& element : element_ranges_) {
    result[element.first] = CalculateRadialGridforAtom(type, element);
  }
  return result;
}

GridContainers::radial_grid EulerMaclaurinGrid::CalculateRadialGridforAtom(
    const std::string& type, const std::pair<std::string, min_exp>& element) {
  GridContainers::radial_grid result;
  Index np = getGridParameters(element.first, type);
  double cutoff = element.second.range;
  result.radius = Eigen::VectorXd::Zero(np);
  result.weight = Eigen::VectorXd::Zero(np);
  double alpha =
      -cutoff /
      (log(1.0 - std::pow((1.0 + double(np)) / (2.0 + double(np)), 3)));
  double factor = 3.0 / (1.0 + double(np));

  for (Index i = 0; i < np; i++) {
    double q = double(i + 1) / (double(np) + 1.0);
    double r = -alpha * std::log(1.0 - std::pow(q, 3));
    double w = factor * alpha * r * r / (1.0 - std::pow(q, 3)) * std::pow(q, 2);
    result.radius[i] = r;
    result.weight[i] = w;
  }
  return result;
}

double EulerMaclaurinGrid::DetermineCutoff(double alpha, Index l, double eps) {
  // determine norm of function
  /* For a function f(r) = r^k*exp(-alpha*r^2) determine
     the radial distance r such that the fraction of the
     function norm that is neglected if the 3D volume
     integration is terminated at a distance r is less
     than or equal to eps. */

  double cutoff = 1.0;     // initial value
  double increment = 0.5;  // increment

  while (increment > 0.01) {
    double residual = CalcResidual(alpha, l, cutoff);
    if (residual > eps) {
      cutoff += increment;
    } else {
      cutoff -= increment;
      if (cutoff < 0.0) {
        cutoff = 0.0;
      }
      increment = 0.5 * increment;
      cutoff += increment;
    }
  }
  return cutoff;
}

double EulerMaclaurinGrid::CalcResidual(double alpha, Index l, double cutoff) {
  return RadialIntegral(alpha, l + 2, cutoff) /
         RadialIntegral(alpha, l + 2, 0.0);
}

double EulerMaclaurinGrid::RadialIntegral(double alpha, Index l,
                                          double cutoff) {
  const double pi = boost::math::constants::pi<double>();
  Index ilo = l % 2;
  double value = 0.0;
  double valexp;
  if (ilo == 0) {
    double expo = std::sqrt(alpha) * cutoff;
    if (expo <= 40.0) {
      value = 0.5 * std::sqrt(pi / alpha) * std::erfc(expo);
    }
  }
  double exponent = alpha * cutoff * cutoff;
  if (exponent > 500.0) {
    valexp = 0.0;
    value = 0.0;
  } else {
    valexp = std::exp(-exponent);
    value = valexp / 2.0 / alpha;
  }
  for (Index i = ilo + 2; i <= l; i += 2) {
    value = (double(i - 1) * value + std::pow(cutoff, i - 1) * valexp) / 2.0 /
            alpha;
  }
  return value;
}

Index EulerMaclaurinGrid::getGridParameters(const std::string& element,
                                            const std::string& type) {
  if (type == "medium") {
    return MediumGrid.at(element);
  } else if (type == "coarse") {
    return CoarseGrid.at(element);
  } else if (type == "xcoarse") {
    return XcoarseGrid.at(element);
  } else if (type == "fine") {
    return FineGrid.at(element);
  } else if (type == "xfine") {
    return XfineGrid.at(element);
  }
  throw std::runtime_error("Grid type " + type + " is not implemented");
  return -1;
}
}  // namespace xtp
}  // namespace votca
