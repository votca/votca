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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE environmentscreening_test

// Standard includes
#include <cmath>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/dipoledipoleinteraction.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_base.h"
#include "votca/xtp/sigmafactory.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vxc_grid.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(environmentscreening_test)

namespace {

QMMolecule Methane() {
  QMMolecule mol("methane", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  return mol;
}

AOBasis MakeBasis(const std::string& name, const QMMolecule& mol) {
  BasisSet bs;
  bs.Load(name);
  AOBasis basis;
  basis.Fill(bs, mol);
  return basis;
}

// All functions of `basis` evaluated at `r`, in AO order.
Eigen::VectorXd AOValuesAt(const AOBasis& basis, const Eigen::Vector3d& r) {
  Eigen::VectorXd v = Eigen::VectorXd::Zero(basis.AOBasisSize());
  for (const AOShell& shell : basis) {
    const AOShell::AOValues vals = shell.EvalAOspace(r);
    v.segment(shell.getStartIndex(), shell.getNumFunc()) = vals.values;
  }
  return v;
}

// Field at each site produced by each auxiliary function taken as a charge
// density, by direct quadrature on a DFT grid -- an independent route to
// the same quantity as AuxFieldAtPoints, sharing nothing with it but the
// basis definition. One pass over the grid: every function is evaluated
// once per point.
Eigen::MatrixXd AuxFieldByQuadrature(
    const Vxc_Grid& grid, const AOBasis& aux,
    const std::vector<Eigen::Vector3d>& sites) {
  Eigen::MatrixXd F =
      Eigen::MatrixXd::Zero(aux.AOBasisSize(), 3 * Index(sites.size()));
  for (Index b = 0; b < grid.getBoxesSize(); ++b) {
    const auto& pts = grid[b].getGridPoints();
    const auto& wts = grid[b].getGridWeights();
    for (std::size_t g = 0; g < pts.size(); ++g) {
      const Eigen::VectorXd chi = wts[g] * AOValuesAt(aux, pts[g]);
      for (std::size_t s = 0; s < sites.size(); ++s) {
        const Eigen::Vector3d d = sites[s] - pts[g];
        const double r = d.norm();
        F.middleCols<3>(3 * Index(s)) += chi * (d / (r * r * r)).transpose();
      }
    }
  }
  return F;
}

// Field at each site of the orbital density psi_n^2, by the same kind of
// quadrature, on a grid built for the DFT basis.
Eigen::VectorXd OrbitalDensityField(const Vxc_Grid& grid, const AOBasis& dft,
                                    const Eigen::VectorXd& coeffs,
                                    const std::vector<Eigen::Vector3d>& sites) {
  Eigen::VectorXd E = Eigen::VectorXd::Zero(3 * Index(sites.size()));
  for (Index b = 0; b < grid.getBoxesSize(); ++b) {
    const auto& pts = grid[b].getGridPoints();
    const auto& wts = grid[b].getGridWeights();
    for (std::size_t g = 0; g < pts.size(); ++g) {
      const double psi = AOValuesAt(dft, pts[g]).dot(coeffs);
      const double w = wts[g] * psi * psi;
      for (std::size_t s = 0; s < sites.size(); ++s) {
        const Eigen::Vector3d d = sites[s] - pts[g];
        const double r = d.norm();
        E.segment<3>(3 * Index(s)) += w * d / (r * r * r);
      }
    }
  }
  return E;
}

PolarSegment SegmentAt(Index id, const std::vector<Eigen::Vector3d>& pos,
                       double alpha) {
  PolarSegment seg("env", id);
  Index k = 0;
  for (const Eigen::Vector3d& p : pos) {
    PolarSite site(k++, "C", p);
    site.setpolarization(alpha * Eigen::Matrix3d::Identity());
    seg.push_back(site);
  }
  return seg;
}

std::vector<Eigen::Vector3d> PositionsOf(
    const std::vector<PolarSegment>& segs) {
  std::vector<Eigen::Vector3d> out;
  for (const PolarSegment& s : segs) {
    for (const PolarSite& site : s) {
      out.push_back(site.getPos());
    }
  }
  return out;
}

// B column by column through the operator PolarRegion itself solves with
// -- DipoleDipoleInteraction under Eigen's CG -- rather than through the
// dense assembly in ReactionFieldKernel, so the two share only
// eeInteractor::FillTholeInteraction.
Eigen::MatrixXd KernelByCG(const Eigen::MatrixXd& F,
                           const std::vector<PolarSegment>& segs,
                           double exp_damp) {
  const eeInteractor interactor(exp_damp);
  DipoleDipoleInteraction A(interactor, segs);
  Eigen::ConjugateGradient<DipoleDipoleInteraction, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setTolerance(1e-14);
  cg.setMaxIterations(10000);
  cg.compute(A);
  Eigen::MatrixXd B(F.rows(), F.rows());
  for (Index P = 0; P < F.rows(); ++P) {
    const Eigen::VectorXd mu = cg.solve(Eigen::VectorXd(F.row(P).transpose()));
    B.col(P) = -F * mu;
  }
  return B;
}

}  // namespace

// Step 1. F against a completely independent evaluation: the auxiliary
// functions sampled on a DFT quadrature grid and the Coulomb field summed
// point by point. Checks the values, the sign convention (field, not
// potential gradient), and that the rows line up with the auxiliary basis
// in the order the three-centre integrals use -- for s, p, d and f
// functions alike.
BOOST_AUTO_TEST_CASE(aux_field_matches_direct_quadrature) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);

  // Generic directions, 11 to 16 bohr out. Far enough that even the most
  // diffuse auxiliary function (exponent 0.207, on carbon) has decayed
  // below 1e-9 at every site: the quadrature reference integrates a
  // 1/r^2 singularity at the site that it does not resolve, and would
  // otherwise be the less accurate of the two. F itself was checked
  // against the closed form down to 2 bohr (see AuxFieldAtPoints).
  const std::vector<Eigen::Vector3d> sites = {Eigen::Vector3d(11.2, 2.1, -1.3),
                                              Eigen::Vector3d(-3.8, 13.6, 5.9),
                                              Eigen::Vector3d(6.3, -8.7, 11.5)};

  const Eigen::MatrixXd F = EnvironmentScreening::AuxFieldAtPoints(aux, sites);
  BOOST_REQUIRE_EQUAL(F.rows(), aux.AOBasisSize());
  BOOST_REQUIRE_EQUAL(F.cols(), 3 * Index(sites.size()));

  Vxc_Grid grid;
  grid.GridSetup("xfine", mol, aux);

  const Eigen::MatrixXd Fq = AuxFieldByQuadrature(grid, aux, sites);

  const double scale = Fq.cwiseAbs().maxCoeff();
  const double err = (F - Fq).cwiseAbs().maxCoeff();
  BOOST_TEST_MESSAGE("max |F|=" << scale << "  max |F - quadrature|=" << err);
  BOOST_CHECK_SMALL(err / scale, 1e-6);
  libint2::finalize();
}

// Step 2. One isolated site has no one to couple to, so its response is
// just alpha: B = -alpha F F^T.
BOOST_AUTO_TEST_CASE(single_site_kernel_is_minus_alpha_F_FT) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  const double alpha = 9.5;
  const std::vector<PolarSegment> segs = {
      SegmentAt(0, {Eigen::Vector3d(6.0, -1.0, 2.5)}, alpha)};
  const Eigen::MatrixXd F =
      EnvironmentScreening::AuxFieldAtPoints(aux, PositionsOf(segs));

  const Eigen::MatrixXd B =
      EnvironmentScreening::ReactionFieldKernel(F, segs, 0.39);
  const Eigen::MatrixXd expected = -alpha * F * F.transpose();
  BOOST_CHECK_SMALL(
      (B - expected).cwiseAbs().maxCoeff() / expected.cwiseAbs().maxCoeff(),
      1e-12);
  libint2::finalize();
}

// Step 2. Coupled sites: the dense factorization against CG on
// DipoleDipoleInteraction, the operator PolarRegion solves. Close enough
// together that the Thole coupling changes the answer substantially --
// checked explicitly, so the comparison cannot pass because the coupling
// happened to be negligible.
BOOST_AUTO_TEST_CASE(coupled_kernel_matches_polarregion_operator) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  const double alpha = 9.5;
  const double exp_damp = 0.39;
  const std::vector<PolarSegment> segs = {
      SegmentAt(0,
                {Eigen::Vector3d(6.5, 0.0, 0.0), Eigen::Vector3d(9.0, 0.3, 0.0),
                 Eigen::Vector3d(7.7, 2.4, 0.2)},
                alpha),
      SegmentAt(
          1, {Eigen::Vector3d(-1.0, 7.0, 2.0), Eigen::Vector3d(-1.2, 9.6, 2.3)},
          alpha)};
  const Eigen::MatrixXd F =
      EnvironmentScreening::AuxFieldAtPoints(aux, PositionsOf(segs));

  const Eigen::MatrixXd B =
      EnvironmentScreening::ReactionFieldKernel(F, segs, exp_damp);
  const Eigen::MatrixXd Bcg = KernelByCG(F, segs, exp_damp);
  const double scale = Bcg.cwiseAbs().maxCoeff();
  BOOST_CHECK_SMALL((B - Bcg).cwiseAbs().maxCoeff() / scale, 1e-9);

  // Symmetric, and negative semidefinite: it screens.
  BOOST_CHECK_SMALL((B - B.transpose()).cwiseAbs().maxCoeff() / scale, 1e-14);
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      B, Eigen::EigenvaluesOnly);
  BOOST_CHECK_LE(es.eigenvalues().maxCoeff(), 1e-12 * scale);

  // The coupling matters here, so the agreement above is not trivial.
  const Eigen::MatrixXd uncoupled =
      EnvironmentScreening::ShellKernel(F, segs, 1.0);
  BOOST_CHECK_GT((B - uncoupled).cwiseAbs().maxCoeff() / scale, 1e-2);
  libint2::finalize();
}

// Step 2. The shell term is the uncoupled response scaled by 1/epsilon.
// At epsilon = 1 it must coincide with the explicit kernel for sites too
// far apart to couple, which ties the two kernels to the same convention.
BOOST_AUTO_TEST_CASE(shell_kernel_is_uncoupled_response_over_epsilon) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  const double alpha = 9.5;
  // Hundreds of bohr apart from one another: coupling ~ alpha/r^3 ~ 1e-6.
  const std::vector<PolarSegment> segs = {
      SegmentAt(0, {Eigen::Vector3d(8.0, 0.0, 0.0)}, alpha),
      SegmentAt(1, {Eigen::Vector3d(0.0, -300.0, 0.0)}, alpha)};
  const Eigen::MatrixXd F =
      EnvironmentScreening::AuxFieldAtPoints(aux, PositionsOf(segs));

  const Eigen::MatrixXd B1 = EnvironmentScreening::ShellKernel(F, segs, 1.0);
  const Eigen::MatrixXd Bx =
      EnvironmentScreening::ReactionFieldKernel(F, segs, 0.39);
  const double scale = Bx.cwiseAbs().maxCoeff();
  BOOST_CHECK_SMALL((B1 - Bx).cwiseAbs().maxCoeff() / scale, 1e-6);

  const Eigen::MatrixXd B4 = EnvironmentScreening::ShellKernel(F, segs, 4.0);
  BOOST_CHECK_SMALL((B4 - 0.25 * B1).cwiseAbs().maxCoeff() / scale, 1e-14);

  BOOST_CHECK_THROW(EnvironmentScreening::ShellKernel(F, segs, 0.5),
                    std::runtime_error);
  libint2::finalize();
}

// Steps 0-3 end to end, with physical meaning. The reaction-field energy
// of an orbital density, 1/2 (nn|v_reac|nn), two ways:
//
//  RI:     1/2 M_nn R M_nn^T, with M from TCMatrix_gwbse and R = T^T B T
//          -- exactly how the GW and BSE code will see it;
//  direct: the field of psi_n^2 at the sites by quadrature, then
//          -1/2 E^T A^-1 E through PolarRegion's own operator.
//
// They share no code path through the aux basis, the metric or the
// kernel. What separates them is the RI fitting error of the density's
// field at distant sites. The energy is monopole-dominated, and the
// Coulomb-metric fit does not conserve charge, so the discrepancy tracks
// the fitted charge of psi_n^2: measured, with aux-def2-svp the HOMO's
// fitted density carries +0.39% charge and the energy is off by 0.60%;
// with aux-def2-qzvp +5.7e-4 and 7.5e-4 (other levels ~1e-5 / ~4e-5).
// That is a property of the auxiliary basis, not of this code. This test
// uses aux-aug-cc-pvtz, all levels within 4.6e-4: nearly as good, and no
// functions above g. qzvp has h functions, which libint2 builds with the
// common --with-max-am=4 refuse (Engine::lmax_exceeded). The exact algebra
// of steps 0-3 is checked separately, to round-off, against the fit
// coefficients themselves.
BOOST_AUTO_TEST_CASE(ri_reaction_energy_matches_direct_evaluation) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis dft = MakeBasis(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/3-21G.xml", mol);
  const AOBasis aux = MakeBasis("aux-aug-cc-pvtz", mol);
  const Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");

  TCMatrix_gwbse tc;
  tc.Initialize(aux.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aux, dft, MOs);

  const double alpha = 9.5;
  const double exp_damp = 0.39;
  // 9.5 bohr and further from the centre, where psi_n^2 has decayed below
  // ~1e-7 of its peak: non-overlapping, as the method assumes, and so the
  // quadrature reference is not integrating across its own singularity.
  // Neighbouring sites 2.5 bohr apart, so the Thole coupling is strong.
  const std::vector<PolarSegment> segs = {
      SegmentAt(
          0,
          {Eigen::Vector3d(9.5, 0.0, 0.0), Eigen::Vector3d(12.0, 0.3, 0.0),
           Eigen::Vector3d(10.7, 2.4, 0.2)},
          alpha),
      SegmentAt(
          1,
          {Eigen::Vector3d(-1.5, 10.0, 2.8), Eigen::Vector3d(-1.7, 12.6, 3.1)},
          alpha)};
  const std::vector<Eigen::Vector3d> sites = PositionsOf(segs);

  const Eigen::MatrixXd F = EnvironmentScreening::AuxFieldAtPoints(aux, sites);
  const Eigen::MatrixXd B =
      EnvironmentScreening::ReactionFieldKernel(F, segs, exp_damp);
  const Eigen::MatrixXd R =
      EnvironmentScreening::SymmetrizedReactionField(B, tc.InvSqrt());

  // Step 3's own guarantees: symmetric, and 1 + R positive definite.
  BOOST_CHECK_SMALL((R - R.transpose()).cwiseAbs().maxCoeff(), 1e-14);
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esR(
      R, Eigen::EigenvaluesOnly);
  // Negative semidefinite up to round-off, which grows with the condition
  // of the metric T (poor for large aux bases), hence relative to the
  // spectrum.
  BOOST_CHECK_LE(esR.eigenvalues().maxCoeff(),
                 1e-9 * std::abs(esR.eigenvalues().minCoeff()));
  BOOST_CHECK_GT(esR.eigenvalues().minCoeff(), -1.0);

  Vxc_Grid grid;
  grid.GridSetup("xfine", mol, dft);
  const eeInteractor interactor(exp_damp);
  DipoleDipoleInteraction A(interactor, segs);
  Eigen::ConjugateGradient<DipoleDipoleInteraction, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setTolerance(1e-14);
  cg.setMaxIterations(10000);
  cg.compute(A);

  for (Index n = 0; n < 5; ++n) {  // every occupied level
    const Eigen::RowVectorXd Mnn = tc[n].row(n);
    const double e_ri = 0.5 * double(Mnn * R * Mnn.transpose());

    // Exact algebra: with fit coefficients c = T M^T (so V^-1 = T T^T),
    // 1/2 M R M^T = 1/2 c^T B c = -1/2 (F^T c)^T A^-1 (F^T c). The field
    // of the fitted density, through PolarRegion's operator; no fit error.
    const Eigen::VectorXd c = tc.InvSqrt() * Mnn.transpose();
    const Eigen::VectorXd Ec = F.transpose() * c;
    const double e_fit = -0.5 * Ec.dot(cg.solve(Ec));
    BOOST_CHECK_SMALL(std::abs(e_ri - e_fit) / std::abs(e_fit), 1e-10);

    const Eigen::VectorXd E = OrbitalDensityField(grid, dft, MOs.col(n), sites);
    const Eigen::VectorXd mu = cg.solve(E);
    const double e_direct = -0.5 * E.dot(mu);

    BOOST_TEST_MESSAGE(
        "level " << n << ": RI " << e_ri << "  direct " << e_direct << "  rel "
                 << std::abs(e_ri - e_direct) / std::abs(e_direct));
    BOOST_CHECK_LT(e_direct, 0.0);  // polarization stabilizes
    BOOST_CHECK_SMALL(std::abs(e_ri - e_direct) / std::abs(e_direct), 1e-3);
  }
  libint2::finalize();
}

// Step 4, physics. The simplest environment there is: one that answers
// any change of charge dQ in the molecule with a uniform potential c dQ,
// v_reac(r, r') = c. Every pair density then feels c times its charge,
// (nm|v_reac|n'm') = c <n|m> <n'|m'> = c delta_nm delta_n'm', and the
// COH+SEX self-energy is exactly
//
//   Sigma^reac_nn' = 1/2 c s_n delta_nn'.
//
// With c < 0, as for any screening environment: occupied levels up by
// |c|/2, virtual levels down by |c|/2 -- the polarization energy of the
// added hole or electron, P = -1/2 (nn|v_reac|nn). SEX alone would give
// the HOMO 2P and the LUMO nothing; this checks the pair is right, and
// that the occupation split falls between HOMO and LUMO.
//
// Built through the real pipeline: aux charges q_P by quadrature, B = c q
// q^T, R = T^T B T. The residual is the RI fit's charge error on the pair
// densities, as in the energy test above.
BOOST_AUTO_TEST_CASE(uniform_reaction_field_shifts_levels_by_plus_minus_P) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis dft = MakeBasis(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/3-21G.xml", mol);
  const AOBasis aux = MakeBasis("aux-aug-cc-pvtz", mol);
  const Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");
  const Index nmo = MOs.cols();
  const Index homo = 4;

  TCMatrix_gwbse tc;
  tc.Initialize(aux.AOBasisSize(), 0, nmo - 1, 0, nmo - 1);
  tc.Fill(aux, dft, MOs);

  // Charges of the auxiliary functions. Only s shells carry any: a pure
  // spherical harmonic with l > 0 integrates to zero, so the quadrature
  // only needs the s shells.
  Vxc_Grid grid;
  grid.GridSetup("xfine", mol, dft);
  Eigen::VectorXd q = Eigen::VectorXd::Zero(aux.AOBasisSize());
  for (Index b = 0; b < grid.getBoxesSize(); ++b) {
    const auto& pts = grid[b].getGridPoints();
    const auto& wts = grid[b].getGridWeights();
    for (std::size_t g = 0; g < pts.size(); ++g) {
      for (const AOShell& shell : aux) {
        if (shell.getL() == L::S) {
          q(shell.getStartIndex()) +=
              wts[g] * shell.EvalAOspace(pts[g]).values(0);
        }
      }
    }
  }
  const double c = -0.1;
  const Eigen::MatrixXd B = c * q * q.transpose();
  const Eigen::MatrixXd R = tc.InvSqrt().transpose() * B * tc.InvSqrt();

  Logger log;
  RPA rpa(log, tc);
  std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("ppm", tc, rpa);
  Sigma_base::options opt;
  opt.homo = homo;
  opt.qpmin = 0;
  opt.qpmax = nmo - 1;
  opt.rpamin = 0;
  opt.rpamax = nmo - 1;
  opt.eta = 1e-3;
  opt.order = 0;
  opt.alpha = 0.0;
  sigma->configure(opt);
  const Eigen::MatrixXd S = sigma->CalcReactionFieldMatrix(R);

  Eigen::VectorXd expected(nmo);
  for (Index n = 0; n < nmo; ++n) {
    expected(n) = 0.5 * c * (n <= homo ? -1.0 : 1.0);
  }
  const Eigen::MatrixXd E = expected.asDiagonal();
  BOOST_TEST_MESSAGE("diag: " << S.diagonal().transpose() << "\nmax dev "
                              << (S - E).cwiseAbs().maxCoeff());
  BOOST_CHECK_GT(S(homo, homo), 0.0);          // HOMO up
  BOOST_CHECK_LT(S(homo + 1, homo + 1), 0.0);  // LUMO down
  BOOST_CHECK_SMALL((S - E).cwiseAbs().maxCoeff() / std::abs(0.5 * c), 5e-3);
  libint2::finalize();
}

namespace {
// Methane with a real Thole environment, as in the RI energy test: the
// fitted integrals, R, and invented but properly ordered RPA energies (the
// identities below hold for any).
struct ScreenedMethane {
  QMMolecule mol = Methane();
  AOBasis dft = MakeBasis(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/3-21G.xml", mol);
  AOBasis aux = MakeBasis("aux-aug-cc-pvtz", mol);
  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");
  Index homo = 4;
  Eigen::VectorXd energies;

  ScreenedMethane() {
    energies.resize(MOs.cols());
    energies << -10.6, -0.75, -0.40, -0.39, -0.38, 0.17, 0.23, 0.24, 0.25, 0.76,
        0.77, 0.78, 1.05, 1.13, 1.14, 1.15, 1.73;
  }

  TCMatrix_gwbse Fill() const {
    TCMatrix_gwbse tc;
    tc.Initialize(aux.AOBasisSize(), 0, MOs.cols() - 1, 0, MOs.cols() - 1);
    tc.Fill(aux, dft, MOs);
    return tc;
  }

  Eigen::MatrixXd TholeR(const TCMatrix_gwbse& tc) const {
    const std::vector<PolarSegment> segs = {
        SegmentAt(
            0,
            {Eigen::Vector3d(7.0, 0.0, 0.0), Eigen::Vector3d(9.5, 0.3, 0.0),
             Eigen::Vector3d(8.2, 2.4, 0.2)},
            9.5),
        SegmentAt(
            1,
            {Eigen::Vector3d(-1.5, 7.5, 2.8), Eigen::Vector3d(-1.7, 10.1, 3.1)},
            9.5)};
    const Eigen::MatrixXd F =
        EnvironmentScreening::AuxFieldAtPoints(aux, PositionsOf(segs));
    return EnvironmentScreening::SymmetrizedReactionField(
        EnvironmentScreening::ReactionFieldKernel(F, segs, 0.39), tc.InvSqrt());
  }

  // The environment that answers any charge with a uniform potential c.
  Eigen::MatrixXd UniformR(const TCMatrix_gwbse& tc, double c) const {
    Vxc_Grid grid;
    grid.GridSetup("xfine", mol, dft);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(aux.AOBasisSize());
    for (Index b = 0; b < grid.getBoxesSize(); ++b) {
      const auto& pts = grid[b].getGridPoints();
      const auto& wts = grid[b].getGridWeights();
      for (std::size_t g = 0; g < pts.size(); ++g) {
        for (const AOShell& shell : aux) {
          if (shell.getL() == L::S) {
            q(shell.getStartIndex()) +=
                wts[g] * shell.EvalAOspace(pts[g]).values(0);
          }
        }
      }
    }
    const Eigen::VectorXd t = tc.InvSqrt().transpose() * q;
    return c * t * t.transpose();
  }

  Eigen::VectorXd SigmaC(TCMatrix_gwbse& tc, const Eigen::VectorXd& w) const {
    Logger log;
    RPA rpa(log, tc);
    rpa.configure(homo, 0, MOs.cols() - 1);
    rpa.setRPAInputEnergies(energies);
    std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("ppm", tc, rpa);
    Sigma_base::options opt;
    opt.homo = homo;
    opt.qpmin = 0;
    opt.qpmax = MOs.cols() - 1;
    opt.rpamin = 0;
    opt.rpamax = MOs.cols() - 1;
    opt.eta = 1e-3;
    opt.order = 0;
    opt.alpha = 0.0;
    sigma->configure(opt);
    sigma->PrepareScreening();
    return sigma->CalcCorrelationDiag(w);
  }
};
}  // namespace

// Step 5, algebra. Dressing M with S = (1 + R)^(1/2) makes the RPA see the
// interaction u = 1 + R: its dielectric matrix becomes 1 + S (eps - 1) S,
// and S eps_dressed^-1 S is the fully screened W = [u^-1 + eps - 1]^-1 --
// the environment's reaction field and the molecule's own screening,
// combined at the level of the Dyson equation, not added.
BOOST_AUTO_TEST_CASE(dressed_rpa_gives_environment_screened_w) {
  libint2::initialize();
  const ScreenedMethane sys;
  TCMatrix_gwbse bare = sys.Fill();
  const Eigen::MatrixXd R = sys.TholeR(bare);
  const Eigen::MatrixXd S = EnvironmentScreening::DressingMatrix(R);
  const Index n = R.rows();
  BOOST_CHECK_SMALL(
      (S * S - Eigen::MatrixXd::Identity(n, n) - R).cwiseAbs().maxCoeff(),
      1e-12);

  Logger log;
  RPA rpa_bare(log, bare);
  rpa_bare.configure(sys.homo, 0, sys.MOs.cols() - 1);
  rpa_bare.setRPAInputEnergies(sys.energies);
  const Eigen::MatrixXd eps = rpa_bare.calculate_epsilon_r(0.0);

  TCMatrix_gwbse dressed = sys.Fill();
  dressed.DressAuxIndex(S);
  RPA rpa_dressed(log, dressed);
  rpa_dressed.configure(sys.homo, 0, sys.MOs.cols() - 1);
  rpa_dressed.setRPAInputEnergies(sys.energies);
  const Eigen::MatrixXd eps_d = rpa_dressed.calculate_epsilon_r(0.0);

  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
  BOOST_CHECK_SMALL((eps_d - (I + S * (eps - I) * S)).cwiseAbs().maxCoeff(),
                    1e-10);
  const Eigen::MatrixXd W = S * eps_d.inverse() * S;
  const Eigen::MatrixXd W_ref = ((I + R).inverse() + eps - I).inverse();
  BOOST_CHECK_SMALL((W - W_ref).cwiseAbs().maxCoeff(), 1e-10);

  // Undressing gives back the bare integrals.
  dressed.UndressAuxIndex();
  for (Index m = 0; m < bare.msize(); ++m) {
    BOOST_CHECK_SMALL((dressed[m] - bare[m]).cwiseAbs().maxCoeff(), 1e-12);
  }
  libint2::finalize();
}

// Step 5, physics. A uniform potential does not couple to a neutral
// excitation -- the transition density psi_i psi_a carries no charge -- so
// an environment that only ever answers with one leaves the molecule's
// dynamic correlation alone: Sigma_c at fixed frequencies is unchanged,
// up to the fitted charge of the transition densities. A real Thole
// environment is not uniform, polarizes in response to the molecule's
// transition dipoles, and does change it -- by far more.
BOOST_AUTO_TEST_CASE(uniform_environment_leaves_correlation_alone) {
  libint2::initialize();
  const ScreenedMethane sys;
  const Eigen::VectorXd w = sys.energies;

  TCMatrix_gwbse bare = sys.Fill();
  const Eigen::MatrixXd R_uniform = sys.UniformR(bare, -0.1);
  const Eigen::MatrixXd R_thole = sys.TholeR(bare);
  const Eigen::VectorXd sc_bare = sys.SigmaC(bare, w);

  TCMatrix_gwbse uni = sys.Fill();
  uni.DressAuxIndex(EnvironmentScreening::DressingMatrix(R_uniform));
  const Eigen::VectorXd sc_uni = sys.SigmaC(uni, w);

  TCMatrix_gwbse thole = sys.Fill();
  thole.DressAuxIndex(EnvironmentScreening::DressingMatrix(R_thole));
  const Eigen::VectorXd sc_thole = sys.SigmaC(thole, w);

  // Occupied levels and the LUMO: away from the (invented) spectrum's
  // PPM poles, where any change is amplified.
  const Index k = sys.homo + 2;
  const double d_uni = (sc_uni - sc_bare).head(k).cwiseAbs().maxCoeff();
  const double d_thole = (sc_thole - sc_bare).head(k).cwiseAbs().maxCoeff();
  BOOST_TEST_MESSAGE("max |dSigma_c|, levels 0.." << k - 1 << ": uniform "
                                                  << d_uni << "  Thole "
                                                  << d_thole);
  BOOST_CHECK_LT(d_uni, 2e-2 * d_thole);
  libint2::finalize();
}

// Step 7. The whole-environment kernel GWBSE builds is the explicit Thole
// region's plus the shell's, each on its own sites, with no coupling.
BOOST_AUTO_TEST_CASE(environment_kernel_is_explicit_plus_shell) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-def2-svp", mol);
  ScreeningEnvironment env;
  env.exp_damp = 0.39;
  env.shell_dielectric = 3.0;
  env.explicit_segments = {SegmentAt(
      0, {Eigen::Vector3d(8.0, 0.0, 0.0), Eigen::Vector3d(10.5, 0.4, 0.1)},
      9.5)};
  env.shell_segments = {SegmentAt(1, {Eigen::Vector3d(0.0, 20.0, 1.0)}, 9.5),
                        SegmentAt(2, {Eigen::Vector3d(-19.0, -3.0, 2.0)}, 9.5)};
  BOOST_CHECK(!env.empty());
  BOOST_CHECK(ScreeningEnvironment().empty());

  const Eigen::MatrixXd B = EnvironmentScreening::Kernel(aux, env);
  // With the default site_width: each part on its own sites, smeared
  // with its own widths.
  const Eigen::MatrixXd Fx = EnvironmentScreening::AuxFieldAtPoints(
      aux, PositionsOf(env.explicit_segments),
      EnvironmentScreening::SiteWidths(env.explicit_segments, env.site_width));
  const Eigen::MatrixXd Fs = EnvironmentScreening::AuxFieldAtPoints(
      aux, PositionsOf(env.shell_segments),
      EnvironmentScreening::SiteWidths(env.shell_segments, env.site_width));
  const Eigen::MatrixXd ref = EnvironmentScreening::ReactionFieldKernel(
                                  Fx, env.explicit_segments, env.exp_damp) +
                              EnvironmentScreening::ShellKernel(
                                  Fs, env.shell_segments, env.shell_dielectric);
  BOOST_CHECK_SMALL((B - ref).cwiseAbs().maxCoeff(), 1e-14);
  BOOST_CHECK_SMALL(EnvironmentScreening::Kernel(aux, ScreeningEnvironment())
                        .cwiseAbs()
                        .maxCoeff(),
                    1e-300);
  libint2::finalize();
}

// Check, the up-front test a QM/MM job runs before its loop, must see the
// same R the GW-BSE run will: Metric builds exactly the T that
// TCMatrix_gwbse::Fill folds into the integrals, and Check's lowest
// eigenvalue is that of SymmetrizedReactionField. And a stable
// environment passes.
BOOST_AUTO_TEST_CASE(check_sees_the_reaction_field_gwbse_uses) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis dft = MakeBasis(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/3-21G.xml", mol);
  const AOBasis aux = MakeBasis("aux-aug-cc-pvtz", mol);
  const Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");
  TCMatrix_gwbse tc;
  tc.Initialize(aux.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aux, dft, MOs);

  const Eigen::MatrixXd T = EnvironmentScreening::Metric(aux);
  BOOST_CHECK_EQUAL((T - tc.InvSqrt()).cwiseAbs().maxCoeff(), 0.0);

  ScreeningEnvironment env;
  env.exp_damp = 0.39;
  env.explicit_segments = {
      SegmentAt(
          0, {Eigen::Vector3d(9.5, 0.0, 0.0), Eigen::Vector3d(12.0, 0.3, 0.0)},
          9.5),
      SegmentAt(1, {Eigen::Vector3d(-1.5, 10.0, 2.8)}, 9.5)};
  env.shell_segments = {SegmentAt(2, {Eigen::Vector3d(0.0, -14.0, 1.0)}, 9.5)};
  env.shell_dielectric = 4.0;

  const Eigen::MatrixXd R = EnvironmentScreening::SymmetrizedReactionField(
      EnvironmentScreening::Kernel(aux, env), tc.InvSqrt());
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      R, Eigen::EigenvaluesOnly);

  const ScreeningCheck check = EnvironmentScreening::Check(aux, mol, env, T);
  BOOST_CHECK(check.ok());
  BOOST_CHECK_EQUAL(check.n_unstable, 0);
  BOOST_CHECK_CLOSE(check.lowest, es.eigenvalues()(0), 1e-8);
  BOOST_CHECK(check.lowest < 0.0);
  BOOST_CHECK(check.report.find("mode 0") != std::string::npos);
  libint2::finalize();
}

// An environment that over-screens. One isotropic site, so R has rank 3
// and its nonzero eigenvalues are -alpha times those of G = F^T T T^T F;
// alpha is chosen so the lowest is exactly -2. Check must say not ok,
// count the unstable modes, and put all of mode 0's reaction on that
// site, which it must name.
BOOST_AUTO_TEST_CASE(check_flags_and_locates_an_unstable_environment) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-aug-cc-pvtz", mol);
  const Eigen::MatrixXd T = EnvironmentScreening::Metric(aux);

  const Eigen::Vector3d where = mol[1].getPos() + Eigen::Vector3d(0, 0, 3.0);
  const Eigen::MatrixXd F =
      EnvironmentScreening::AuxFieldAtPoints(aux, {where});
  const Eigen::MatrixXd TF = T.transpose() * F;
  // Dynamic size on purpose: GCC's -Wmaybe-uninitialized misfires on the
  // fixed-size 3x3 solver's eigenvalue storage.
  const Eigen::MatrixXd G = TF.transpose() * TF;
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eg(G);
  const double alpha = 2.0 / eg.eigenvalues().maxCoeff();

  ScreeningEnvironment env;
  env.site_width = 0.0;  // point sites: the construction above is for them
  env.explicit_segments = {SegmentAt(7, {where}, alpha)};
  const ScreeningCheck check = EnvironmentScreening::Check(aux, mol, env, T);

  BOOST_CHECK(!check.ok());
  BOOST_CHECK_CLOSE(check.lowest, -2.0, 1e-6);
  BOOST_CHECK_GE(check.n_unstable, 1);
  BOOST_CHECK(check.report.find("100.0%  segment 7") != std::string::npos);
  BOOST_CHECK(check.report.find("from QM atom 1 H") != std::string::npos);
  BOOST_TEST_MESSAGE(check.report);

  // The same site smeared at the default width no longer over-screens by
  // a factor of two: the mode was the site sitting in the tails of the
  // diffuse functions, which is what the smearing is for.
  env.site_width = 0.5;
  const ScreeningCheck smeared = EnvironmentScreening::Check(aux, mol, env, T);
  BOOST_CHECK(smeared.lowest > check.lowest + 0.5);
  BOOST_TEST_MESSAGE(smeared.report);
  libint2::finalize();
}

// Smeared sites. A spherical charge acts as a point outside its own
// extent, so far from the auxiliary functions a Gaussian site feels the
// field a point feels; and as the width goes to zero the two coincide
// everywhere. Both checked for s, p, d and f functions of a diffuse
// basis, at a close and a distant site.
BOOST_AUTO_TEST_CASE(smeared_site_field_limits) {
  libint2::initialize();
  const QMMolecule mol = Methane();
  const AOBasis aux = MakeBasis("aux-aug-cc-pvtz", mol);
  const Eigen::Vector3d H = mol[1].getPos();
  const std::vector<Eigen::Vector3d> close = {H + Eigen::Vector3d(0, 0, 3.0)};
  const std::vector<Eigen::Vector3d> far = {Eigen::Vector3d(0, 0, 30.0)};

  auto rel = [](const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
    return (a - b).norm() / b.norm();
  };
  const Eigen::MatrixXd Fclose =
      EnvironmentScreening::AuxFieldAtPoints(aux, close);
  const Eigen::MatrixXd Ffar = EnvironmentScreening::AuxFieldAtPoints(aux, far);

  // Far: smearing over 1.5 bohr makes no difference.
  BOOST_CHECK_SMALL(
      rel(EnvironmentScreening::AuxFieldAtPoints(aux, far, {1.5}), Ffar), 1e-8);
  // Width -> 0: point sites recovered, also inside the tails, and at the
  // rate a smeared charge must approach a point: the difference is the
  // Gaussian's second moment (R^2/2 per direction)
  // acting on the local density of chi_Q, so it goes as R^2.
  const double d1 =
      rel(EnvironmentScreening::AuxFieldAtPoints(aux, close, {1e-2}), Fclose);
  const double d2 =
      rel(EnvironmentScreening::AuxFieldAtPoints(aux, close, {5e-3}), Fclose);
  BOOST_CHECK_CLOSE(d1 / d2, 4.0, 1.0);
  BOOST_CHECK_SMALL(d2, 1e-5);
  // Inside the tails a finite width does change F -- that is the point.
  BOOST_CHECK(rel(EnvironmentScreening::AuxFieldAtPoints(aux, close, {1.5}),
                  Fclose) > 1e-3);
  // A width per point, mixed zero and non-zero, is honoured point by point.
  const Eigen::MatrixXd mixed = EnvironmentScreening::AuxFieldAtPoints(
      aux, {close[0], far[0]}, {0.0, 1.5});
  BOOST_CHECK_EQUAL((mixed.leftCols<3>() - Fclose).cwiseAbs().maxCoeff(), 0.0);
  BOOST_CHECK_SMALL(rel(mixed.rightCols<3>(), Ffar), 1e-8);
  BOOST_CHECK_THROW(
      EnvironmentScreening::AuxFieldAtPoints(aux, close, {1.0, 2.0}),
      std::runtime_error);

  // SiteWidths: site_width * (tr(alpha)/3)^(1/3).
  const std::vector<double> w =
      EnvironmentScreening::SiteWidths({SegmentAt(0, {close[0]}, 8.0)}, 0.5);
  BOOST_REQUIRE_EQUAL(w.size(), 1);
  BOOST_CHECK_CLOSE(w[0], 1.0, 1e-10);
  BOOST_CHECK_THROW(
      EnvironmentScreening::SiteWidths({SegmentAt(0, {close[0]}, 8.0)}, -1.0),
      std::runtime_error);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
