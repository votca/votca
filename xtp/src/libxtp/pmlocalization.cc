/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
 * Reference- A fast intrinsic localization procedure applicable for ab initio
 * and semiempirical linear combination of atomic orbital wave functions J.
 * Chem. Phys. 90, 4916 (1989); https://doi.org/10.1063/1.456588 János Pipek and
 * Paul G. Mezey
 */

#include "votca/xtp/pmlocalization.h"
#include "votca/xtp/aomatrix.h"
#include <limits>

namespace votca {
namespace xtp {

void PMLocalization::computePML(Orbitals &orbitals) {

  if (method_ == "jacobi-sweeps") {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Using Jacobi-Sweeps" << std::flush;
    computePML_JS(orbitals);
  } else if (method_ == "unitary-optimizer") {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Using Unitary Optimizer" << std::flush;
    computePML_UT(orbitals);
  }
}

double PMLocalization::cost(const Eigen::MatrixXd &W,
                            const std::vector<Eigen::MatrixXd> &Sat_all,
                            const Index nat) const {

  double Dinv = 0.0;
  double p = 2.0;  // standard PM
  for (Index iat = 0; iat < nat; iat++) {
    Eigen::MatrixXd qw = Sat_all[iat] * W;
    for (Index i = 0; i < W.cols(); i++) {
      double Qa = W.col(i).transpose() * qw.col(i);
      Dinv += std::pow(Qa, p);
    }
  }
  return Dinv;
}

std::pair<double, Eigen::MatrixXd> PMLocalization::cost_derivative(
    const Eigen::MatrixXd &W, const std::vector<Eigen::MatrixXd> &Sat_all,
    const Index nat) const {
  Eigen::MatrixXd Jderiv = Eigen::MatrixXd::Zero(n_occs_, n_occs_);
  double Dinv = 0.0;
  double p = 2.0;  // standard PM
  for (Index iat = 0; iat < nat; iat++) {
    Eigen::MatrixXd qw = Sat_all[iat] * W;
    for (Index i = 0; i < W.cols(); i++) {
      double qwp = W.col(i).transpose() * qw.col(i);
      Dinv += std::pow(qwp, p);
      double t = p * std::pow(qwp, p - 1);
      for (Index j = 0; j < W.cols(); j++) {
        Jderiv(j, i) += t * qw(j, i);
      }
    }
  }
  return {Dinv, Jderiv};
}

Eigen::VectorXd PMLocalization::fit_polynomial(const Eigen::VectorXd &x,
                                               const Eigen::VectorXd &y) const {

  // Fit function to polynomial

  if (x.size() != y.size()) {
    throw std::runtime_error("x and y have different dimensions!\n");
  }
  Index N = x.size();
  Index deg = N;

  // Form matrix
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, deg);
  for (Index i = 0; i < N; i++) {
    for (Index j = 0; j < deg; j++) {
      A(i, j) = std::pow(x(i), j);
    }
  }
  // Solve for coefficients: A * c = y
  return A.colPivHouseholderQr().solve(y);
}

double PMLocalization::find_smallest_step(const Eigen::VectorXd &coeff) const {

  // get the complex roots of the polynomial
  Eigen::VectorXcd complex_roots =
      find_complex_roots(coeff.cast<std::complex<double>>());

  // Real roots
  std::vector<double> real_roots;
  for (Index i = 0; i < complex_roots.size(); i++) {
    if (std::abs(std::imag(complex_roots(i))) <
        10 * std::numeric_limits<double>::epsilon()) {
      real_roots.push_back(std::real(complex_roots(i)));
    }
  }

  // Sort roots
  std::sort(real_roots.begin(), real_roots.end());

  double step = 0.0;
  // for (Index i = 0; i < Index(real_roots.size()); i++) {
  for (auto root : real_roots) {
    // Omit extremely small steps because they might get you stuck.
    if (root > std::sqrt(std::numeric_limits<double>::epsilon())) {
      step = root;
      break;
    }
  }
  return step;
}

Eigen::VectorXcd PMLocalization::find_complex_roots(
    const Eigen::VectorXcd &coeff) const {

  // Coefficient of highest order term must be nonzero.
  Index order = coeff.size();
  while (coeff(order - 1) == 0.0 && order >= 1) order--;

  if (order == 1) {
    // Zeroth degree - no zeros!
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Polynomial is constant - no zero can be found."
        << std::flush;
  }

  // Form companion matrix
  Eigen::MatrixXcd cmat = companion_matrix(coeff.head(order));

  // and diagonalize it
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(cmat);

  // Return the roots (unsorted)
  return es.eigenvalues();
}

Eigen::MatrixXcd PMLocalization::companion_matrix(
    const Eigen::VectorXcd &c) const {
  if (c.size() <= 1) {
    // Dummy return
    Eigen::MatrixXcd dum;
    return dum;
  }

  // Form companion matrix
  Index N = c.size() - 1;
  if (c(N) == 0.0) {
    throw std::runtime_error("Coefficient of highest term vanishes!\n");
  }

  Eigen::MatrixXcd companion = Eigen::MatrixXcd::Zero(N, N);

  // First row - coefficients normalized to that of highest term.
  for (Index j = 0; j < N; j++) {
    companion(0, j) = -c(N - (j + 1)) / c(N);
  }

  // Fill out the unit matrix part
  for (Index j = 1; j < N; j++) {
    companion(j, j - 1) = 1.0;
  }

  return companion;
}

Eigen::MatrixXd PMLocalization::rotate_W(const double step,
                                         const Eigen::MatrixXd &W,
                                         const Eigen::VectorXcd &eval,
                                         const Eigen::MatrixXcd &evec) const {

  Eigen::VectorXcd temp = (step * eval).array().exp();
  Eigen::MatrixXd W_rotated =
      (evec * temp.asDiagonal() * evec.adjoint()).real() * W;  //.real because
                                                               // we assume
                                                               // only real W

  return W_rotated;
}

void PMLocalization::computePML_UT(Orbitals &orbitals) {

  Eigen::MatrixXd occupied_orbitals = orbitals.MOs().eigenvectors().leftCols(
      orbitals.getNumberOfAlphaElectrons());

  // initialize a unitary matrix (as Idenity matrix for now)
  n_occs_ = orbitals.getNumberOfAlphaElectrons();
  W_ = Eigen::MatrixXd::Identity(n_occs_, n_occs_);

  // prepare Mulliken charges
  // get overlap matrix
  aobasis_ = orbitals.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis_);
  overlap_ = overlap.Matrix();
  numfuncpatom_ = aobasis_.getFuncPerAtom();
  Index numatoms = Index(numfuncpatom_.size());

  // could be a bit memory expensive
  std::vector<Eigen::MatrixXd> Sat_all = setup_pop_matrices(occupied_orbitals);
  XTP_LOG(Log::info, log_) << TimeStamp() << " Calculated charge matrices"
                           << std::flush;

  // initialize Riemannian gradient and search direction matrices
  G_ = Eigen::MatrixXd::Zero(n_occs_, n_occs_);
  H_ = Eigen::MatrixXd::Zero(n_occs_, n_occs_);

  bool converged = false;
  Index iteration = 0;
  std::string update_type;
  while (!converged) {
    // Store old gradient and search direction
    G_old_ = G_;
    H_old_ = H_;

    // calculate cost and its derivative wrt unitary matrix for current W
    auto [J, Jderiv] = cost_derivative(W_, Sat_all, numatoms);
    J_ = J;
    XTP_LOG(Log::info, log_)
        << TimeStamp() << " Calculated cost function and its W-derivative"
        << std::flush;

    // calculate Riemannian derivative
    G_ = Jderiv * W_.transpose() - W_ * Jderiv.transpose();

    if (iteration == 0 || (iteration - 1) % W_.cols() == 0) {
      // calculate search direction using SDSA for now
      update_type = "SDSA";
      H_ = G_;
    } else {
      // calculate search direction from CGPR update
      update_type = "CGPR";
      double gamma = inner_prod(G_ - G_old_, G_) / inner_prod(G_old_, G_old_);
      H_ = G_ + gamma * H_old_;
      // make sure H_ is exactly skew symmetric
      H_ = 0.5 * (H_ - H_.transpose());

      // Check that update is OK
      if (inner_prod(G_, H_) < 0.0) {
        H_ = G_;
        XTP_LOG(Log::error, log_)
            << TimeStamp() << "CG search direction reset." << std::flush;
      }
    }

    Index orderW = 4;  // for PM
    // H is skew symmetric, real so should have purely imaginary eigenvalues
    // in pairs +-eval, and 0, if dim is odd.
    Eigen::EigenSolver<Eigen::MatrixXd> es(H_);
    Eigen::VectorXcd Hval = es.eigenvalues();
    Eigen::MatrixXcd Hvec = es.eigenvectors();
    double wmax = Hval.cwiseAbs().maxCoeff();
    double Tmu = 2.0 * tools::conv::Pi / (static_cast<double>(orderW) * wmax);
    double step;
    // line optimization via polynomial fit
    Index npoints = 4;  // to be adapted/adaptable
    double deltaTmu = Tmu / static_cast<double>(npoints - 1);
    int halved = 0;

    // finding optimal step
    while (true) {

      Eigen::VectorXd mu(npoints);
      Eigen::VectorXd derivative_points(npoints);  // cost function derivative
      Eigen::VectorXd cost_points(npoints);        // cost function value

      for (Index i = 0; i < npoints; i++) {
        mu(i) = static_cast<double>(i) * deltaTmu;
        Eigen::MatrixXd W_rotated = rotate_W(mu(i), W_, Hval, Hvec);

        // calculate cost and derivative for this rotated W matrix
        auto [cost, der] = cost_derivative(W_rotated, Sat_all, numatoms);
        cost_points(i) = cost;
        derivative_points(i) =
            2.0 *
            std::real((der * W_rotated.transpose() * H_.transpose()).trace());
      }

      // Check sign of the derivative
      if (derivative_points(0) < 0.0) {
        XTP_LOG(Log::error, log_)
            << TimeStamp() << "Derivative is negative!" << mu << "  "
            << derivative_points << std::flush;
      }

      // Fit to polynomial of order p
      Eigen::VectorXd polyfit_coeff = fit_polynomial(mu, derivative_points);

      // Find step as smallest real zero of the polynomial
      step = find_smallest_step(polyfit_coeff);

      // is step too far?
      if (step > 0.0 && step <= Tmu) {
        // is in range, let's continue
        Eigen::MatrixXd W_new = rotate_W(step, W_, Hval, Hvec);

        // has objective function value changed in the right direction?
        double J_new = cost(W_new, Sat_all, numatoms);
        double delta_J = J_new - J_;

        if (delta_J > 0.0) {
          //  we accept and update
          W_old_ = W_;
          W_ = W_new;
          J_old_ = J_;
          J_ = J_new;
          break;
        } else {
          XTP_LOG(Log::error, log_)
              << TimeStamp()
              << "    WARNING: Cost function is decreasing. deltaJ = "
              << delta_J << std::flush;

          if (halved < 10) {
            XTP_LOG(Log::error, log_)
                << TimeStamp() << "Trying halved maximum step size "
                << std::flush;
            halved++;
            deltaTmu /= 2.0;
            continue;
          } else {
            throw std::runtime_error(
                "Problem in polynomial line search - could not find suitable "
                "extremum!\n");
          }
        }
      } else {
        // now adjust step search, if original step went too far
        XTP_LOG(Log::error, log_)
            << TimeStamp()
            << "Step went beyond max step, trying reduced max step..."
            << std::flush;
        if (halved < 4) {
          halved++;
          deltaTmu /= 2.0;
          continue;
        } else {
          throw std::runtime_error(
              "Problem in polynomial line search - could not find suitable "
              "extremum!\n");
        }
      }
    }

    double G_norm = inner_prod(G_, G_);

    XTP_LOG(Log::error, log_)
        << (boost::format(" UT iteration = %1$6i (%6$4.s) Tmu = %4$4.2e mu_opt "
                          "= %5$1.4f |deltaJ| = %2$4.2e |G| = %3$4.2e") %
            (iteration) % std::abs(J_ - J_old_) % G_norm % Tmu % step %
            update_type)
               .str()
        << std::flush;

    if (iteration > 0 &&
        (G_norm < G_threshold_ && std::abs(J_ - J_old_) < J_threshold_)) {
      converged = true;
    }
    iteration++;
  }

  // all done, what are the actual LMOS?
  localized_orbitals_ =
      (W_.transpose() * occupied_orbitals.transpose()).transpose();  //?
  check_orthonormality();
  // determine the energies of the localized orbitals
  Eigen::VectorXd energies = calculate_lmo_energies(orbitals);

  // sort the LMOS according to energies
  auto [LMOS, LMOS_energies] = sort_lmos(energies);
  // print info
  for (Index i = 0; i < LMOS_energies.size(); i++) {
    XTP_LOG(Log::error, log_)
        << (boost::format(" LMO index = %1$4i Energy = %2$10.5f Hartree ") %
            (i) % LMOS_energies(i))
               .str()
        << std::flush;
  }
  orbitals.setLMOs(LMOS);
  orbitals.setLMOs_energies(LMOS_energies);
}

void PMLocalization::computePML_JS(Orbitals &orbitals) {
  // initialize with occupied CMOs
  localized_orbitals_ = orbitals.MOs().eigenvectors().leftCols(
      orbitals.getNumberOfAlphaElectrons());
  aobasis_ = orbitals.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis_);
  overlap_ = overlap.Matrix();

  XTP_LOG(Log::error, log_) << std::flush;
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Starting localization of orbitals" << std::flush;

  // determine initial penalty_function
  initial_penalty();

  Index iteration = 1;

  while (iteration < nrOfIterations_) {

    XTP_LOG(Log::info, log_) << "Iteration: " << iteration << std::flush;

    Index maxrow, maxcol;
    double max_penalty = PM_penalty_.maxCoeff(&maxrow, &maxcol);

    XTP_LOG(Log::info, log_)
        << "maximum of penalty function: " << max_penalty << std::flush;

    if (max_penalty < convergence_limit_) break;

    XTP_LOG(Log::info, log_)
        << "Orbitals to be changed: " << maxrow << " " << maxcol << std::flush;

    Eigen::MatrixX2d max_orbs(localized_orbitals_.rows(), 2);
    max_orbs << localized_orbitals_.col(maxrow),
        localized_orbitals_.col(maxcol);
    Eigen::MatrixX2d rotated_orbs = rotateorbitals(max_orbs, maxrow, maxcol);
    localized_orbitals_.col(maxrow) = rotated_orbs.col(0);
    localized_orbitals_.col(maxcol) = rotated_orbs.col(1);

    update_penalty(maxrow, maxcol);

    iteration++;
  }
  if (iteration == nrOfIterations_) {
    throw std::runtime_error(
        "Localization with Jacobi-Sweeps did not converge");
  }
  XTP_LOG(Log::error, log_) << TimeStamp() << " Orbitals localized after "
                            << iteration + 1 << " iterations" << std::flush;

  // check if localized orbitals orthonormal, if nor warn
  check_orthonormality();

  // determine the energies of the localized orbitals
  Eigen::VectorXd energies = calculate_lmo_energies(orbitals);

  // sort the LMOS according to energies
  auto [LMOS, LMOS_energies] = sort_lmos(energies);
  // print info
  for (Index i = 0; i < LMOS_energies.size(); i++) {
    XTP_LOG(Log::error, log_)
        << (boost::format(" LMO index = %1$4i Energy = %2$10.5f Hartree ") %
            (i) % LMOS_energies(i))
               .str()
        << std::flush;
  }
  orbitals.setLMOs(LMOS);
  orbitals.setLMOs_energies(LMOS_energies);
}

// sort the LMOs according to their energy
std::pair<Eigen::MatrixXd, Eigen::VectorXd> PMLocalization::sort_lmos(
    const Eigen::VectorXd &energies) {
  Eigen::VectorXd LMOS_energies(energies.size());
  Eigen::MatrixXd LMOS(localized_orbitals_.rows(), localized_orbitals_.cols());

  // sort the LMOs according to energy
  std::vector<std::pair<double, int>> vp;

  // Inserting element in pair vector
  // to keep track of previous indexes
  for (int i = 0; i < energies.size(); ++i) {
    vp.push_back(std::make_pair(energies(i), i));
  }
  // Sorting pair vector
  std::sort(vp.begin(), vp.end());

  for (Index i = 0; i < energies.size(); i++) {
    LMOS_energies(i) = vp[i].first;
    LMOS.col(i) = localized_orbitals_.col(vp[i].second);
  }

  return {LMOS, LMOS_energies};
}

// calculate energies of LMOs
Eigen::VectorXd PMLocalization::calculate_lmo_energies(
    const Orbitals &orbitals) {
  Eigen::MatrixXd h = overlap_ * orbitals.MOs().eigenvectors() *
                      orbitals.MOs().eigenvalues().asDiagonal() *
                      orbitals.MOs().eigenvectors().transpose() * overlap_;

  return (localized_orbitals_.transpose() * h * localized_orbitals_).diagonal();
}

// check orthonormality of LMOs
void PMLocalization::check_orthonormality() {

  Eigen::MatrixXd norm =
      localized_orbitals_.transpose() * overlap_ * localized_orbitals_;
  Eigen::MatrixXd check_norm =
      norm - Eigen::MatrixXd::Identity(norm.rows(), norm.cols());
  bool not_orthonormal = (check_norm.cwiseAbs().array() > 1e-5).any();
  if (not_orthonormal) {
    XTP_LOG(Log::error, log_) << TimeStamp()
                              << " WARNING: Localized orbtials are not "
                                 "orthonormal. Proceed with caution! "
                              << std::flush;
    XTP_LOG(Log::info, log_) << TimeStamp() << " LMOs * S * LMOs" << std::flush;
    for (Index i = 0; i < norm.rows(); i++) {
      for (Index j = 0; j < norm.cols(); j++) {
        if (std::abs(check_norm(i, j)) > 1e-5) {
          XTP_LOG(Log::info, log_)
              << TimeStamp()
              << (boost::format("  Element (%1$4i,%2$4i) = %3$8.2e") % (i) % j %
                  norm(i, j))
                     .str()
              << std::flush;
        }
      }
    }
  }

  return;
}

// Function to rotate the 2 maximum orbitals (s and t)
Eigen::MatrixX2d PMLocalization::rotateorbitals(const Eigen::MatrixX2d &maxorbs,
                                                const Index s, const Index t) {
  const double gamma =
      0.25 *
      asin(B_(s, t) / sqrt((A_(s, t) * A_(s, t)) + (B_(s, t) * B_(s, t))));
  Eigen::MatrixX2d rotatedorbitals(maxorbs.rows(), 2);
  rotatedorbitals.col(0) =
      (std::cos(gamma) * maxorbs.col(0)) + (std::sin(gamma) * maxorbs.col(1));
  rotatedorbitals.col(1) = -1 * (std::sin(gamma) * maxorbs.col(0)) +
                           (std::cos(gamma) * maxorbs.col(1));
  XTP_LOG(Log::info, log_) << "Sine of the rotation angle = " << std::sin(gamma)
                           << std::flush;
  return rotatedorbitals;
}

Eigen::VectorXd PMLocalization::pop_per_atom(const Eigen::VectorXd &orbital) {

  Eigen::RowVectorXd MullikenPop_orb_per_basis =
      (orbital.asDiagonal() * overlap_ * orbital.asDiagonal()).colwise().sum();
  Index start = 0;

  Eigen::VectorXd per_atom = Eigen::VectorXd::Zero(Index(numfuncpatom_.size()));
  for (Index atom_id = 0; atom_id < Index(numfuncpatom_.size()); atom_id++) {
    per_atom(atom_id) =
        MullikenPop_orb_per_basis.segment(start, numfuncpatom_[atom_id]).sum();
    start += numfuncpatom_[atom_id];
  }

  return per_atom;
}

// Determine PM cost function based on Mulliken populations
void PMLocalization::initial_penalty() {

  PM_penalty_ = Eigen::MatrixXd::Zero(localized_orbitals_.cols(),
                                      localized_orbitals_.cols());
  // Variable names A and B are used directly as described in the paper above
  A_ = Eigen::MatrixXd::Zero(localized_orbitals_.cols(),
                             localized_orbitals_.cols());
  B_ = A_;

  numfuncpatom_ = aobasis_.getFuncPerAtom();

  // get the s-s elements first ("diagonal in orbital")
  MullikenPop_orb_per_atom_ = Eigen::MatrixXd::Zero(
      localized_orbitals_.cols(), Index(numfuncpatom_.size()));
#pragma omp parallel for
  for (Index s = 0; s < localized_orbitals_.cols(); s++) {
    MullikenPop_orb_per_atom_.row(s) = pop_per_atom(localized_orbitals_.col(s));
  }

// now we only need to calculate the off-diagonals explicitly
#pragma omp parallel for
  for (Index s = 0; s < localized_orbitals_.cols(); s++) {
    Eigen::MatrixXd s_overlap =
        localized_orbitals_.col(s).asDiagonal() * overlap_;

    for (Index t = s + 1; t < localized_orbitals_.cols(); t++) {

      Eigen::Vector2d temp = offdiag_penalty_elements(s_overlap, s, t);

      A_(s, t) = temp(0);
      B_(s, t) = temp(1);
      PM_penalty_(s, t) =
          A_(s, t) + sqrt((A_(s, t) * A_(s, t)) + (B_(s, t) * B_(s, t)));
    }
  }
  return;
}

std::vector<Eigen::MatrixXd> PMLocalization::setup_pop_matrices(
    const Eigen::MatrixXd &occ_orbitals) {

  // initialize everything
  Index numatoms = Index(numfuncpatom_.size());
  Index noccs = occ_orbitals.cols();

  std::vector<Eigen::MatrixXd> Qat;
  for (Index iat = 0; iat < numatoms; iat++) {
    Qat.push_back(Eigen::MatrixXd::Zero(noccs, noccs));
  }

  numfuncpatom_ = aobasis_.getFuncPerAtom();

  // get the s-s elements first ("diagonal in orbital")
  MullikenPop_orb_per_atom_ =
      Eigen::MatrixXd::Zero(noccs, Index(numfuncpatom_.size()));
#pragma omp parallel for
  for (Index s = 0; s < noccs; s++) {
    MullikenPop_orb_per_atom_.row(s) = pop_per_atom(occ_orbitals.col(s));
  }

  // put this on the diagonals
  for (Index iat = 0; iat < numatoms; iat++) {
    Qat[iat].diagonal() = MullikenPop_orb_per_atom_.col(iat);
  }

// now fill the offdiagonals
#pragma omp parallel for
  for (Index s = 0; s < noccs; s++) {

    Eigen::MatrixXd s_overlap = occ_orbitals.col(s).asDiagonal() * overlap_;

    for (Index t = s + 1; t < noccs; t++) {

      Eigen::MatrixXd splitwiseMullikenPop_orb_SandT =
          s_overlap * occ_orbitals.col(t).asDiagonal();
      Eigen::RowVectorXd MullikenPop_orb_SandT_per_basis =
          0.5 * (splitwiseMullikenPop_orb_SandT.colwise().sum() +
                 splitwiseMullikenPop_orb_SandT.rowwise().sum().transpose());

      Index start = 0;

      for (Index atom_id = 0; atom_id < Index(numfuncpatom_.size());
           atom_id++) {
        Qat[atom_id](s, t) = MullikenPop_orb_SandT_per_basis
                                 .segment(start, numfuncpatom_[atom_id])
                                 .sum();
        Qat[atom_id](t, s) = Qat[atom_id](s, t);

        start += numfuncpatom_[atom_id];
      }
    }
  }

  return Qat;
}

Eigen::Vector2d PMLocalization::offdiag_penalty_elements(
    const Eigen::MatrixXd &s_overlap, Index s, Index t) {

  Eigen::MatrixXd splitwiseMullikenPop_orb_SandT =
      s_overlap * localized_orbitals_.col(t).asDiagonal();
  Eigen::RowVectorXd MullikenPop_orb_SandT_per_basis =
      0.5 * (splitwiseMullikenPop_orb_SandT.colwise().sum() +
             splitwiseMullikenPop_orb_SandT.rowwise().sum().transpose());

  Index start =
      0;  // This helps to sum only over the basis functions on an atom
  double Ast = 0;
  double Bst = 0;

  for (Index atom_id = 0; atom_id < Index(numfuncpatom_.size()); atom_id++) {
    double MullikenPop_orb_SandT_per_atom =
        MullikenPop_orb_SandT_per_basis.segment(start, numfuncpatom_[atom_id])
            .sum();

    Ast += MullikenPop_orb_SandT_per_atom * MullikenPop_orb_SandT_per_atom -
           0.25 * ((MullikenPop_orb_per_atom_(s, atom_id) -
                    MullikenPop_orb_per_atom_(t, atom_id)) *
                   (MullikenPop_orb_per_atom_(s, atom_id) -
                    MullikenPop_orb_per_atom_(t, atom_id)));

    Bst += MullikenPop_orb_SandT_per_atom *
           (MullikenPop_orb_per_atom_(s, atom_id) -
            MullikenPop_orb_per_atom_(t, atom_id));
    start += numfuncpatom_[atom_id];
  }

  Eigen::Vector2d out(Ast, Bst);

  return out;
}

// Update PM cost function based on Mulliken populations after rotations
void PMLocalization::update_penalty(Index orb1, Index orb2) {

  // update the get the s-s elements for orb1 and orb2
#pragma omp parallel for
  for (Index s = 0; s < localized_orbitals_.cols(); s++) {
    if (s == orb1 || s == orb2) {

      MullikenPop_orb_per_atom_.row(s) =
          pop_per_atom(localized_orbitals_.col(s));
    }
  }

// now we only need to calculate the off-diagonals explicitly for all
// pairs involving orb1 or orb2
#pragma omp parallel for
  for (Index s = 0; s < localized_orbitals_.cols(); s++) {
    Eigen::MatrixXd s_overlap =
        localized_orbitals_.col(s).asDiagonal() * overlap_;

    for (Index t = s + 1; t < localized_orbitals_.cols(); t++) {

      // we do this only if any of s or t matches orb1 or orb2
      if (s == orb1 || s == orb2 || t == orb1 || t == orb2) {

        Eigen::Vector2d temp = offdiag_penalty_elements(s_overlap, s, t);
        A_(s, t) = temp(0);
        B_(s, t) = temp(1);
        PM_penalty_(s, t) =
            A_(s, t) + sqrt((A_(s, t) * A_(s, t)) + (B_(s, t) * B_(s, t)));
      }
    }
  }
  return;
}

}  // namespace xtp
}  // namespace votca
