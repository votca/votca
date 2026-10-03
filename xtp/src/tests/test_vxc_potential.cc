/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE vxc_potential_test

// Standard includes
#include <array>
#include <fstream>
#include <sstream>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_functionals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"
#include "xtp_libint2.h"
#include <votca/tools/eigenio_matrixmarket.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(vxc_potential_test)

AOBasis CreateBasis(const QMMolecule& mol) {
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

Eigen::MatrixXd DMat() {
  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/dmat.mm");
  return dmat;
}

BOOST_AUTO_TEST_CASE(vxc_test) {
  libint2::initialize();
  QMMolecule mol("none", 0);

  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/vxc_potential/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = DMat();
  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);
  Vxc_Potential<Vxc_Grid> num(grid);
  num.setXCfunctional("XC_GGA_X_PBE XC_GGA_C_PBE");

  BOOST_CHECK_EQUAL(grid.getGridSize(), grid.getGridpoints().size());
  BOOST_CHECK_EQUAL(grid.getGridSize(), 53404);
  BOOST_CHECK_EQUAL(grid.getBoxesSize(), 51);

  BOOST_CHECK_CLOSE(num.getExactExchange("XC_GGA_X_PBE XC_GGA_C_PBE"), 0.0,
                    1e-5);
  BOOST_CHECK_CLOSE(num.getExactExchange("XC_HYB_GGA_XC_PBEH"), 0.25, 1e-5);
  Mat_p_Energy e_vxc = num.IntegrateVXC(dmat);

  Eigen::MatrixXd vxc_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/vxc_ref.mm");

  bool check_vxc = e_vxc.matrix().isApprox(vxc_ref, 0.0001);

  BOOST_CHECK_CLOSE(e_vxc.energy(), -4.6303432151572643, 1e-5);
  if (!check_vxc) {
    std::cout << "ref" << std::endl;
    std::cout << vxc_ref << std::endl;
    std::cout << "calc" << std::endl;
    std::cout << e_vxc.matrix() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_vxc, 1);

  libint2::finalize();
}

namespace {
// The point-by-point integration the batched code replaced, calling libxc
// directly: one AO evaluation, one libxc call and one rank-1 update per grid
// point. Reference for checking the batched IntegrateVXC/IntegrateVXCSpin to
// rounding.
struct Libxc {
  Libxc(const std::string& functional, int polarization) {
    Vxc_Functionals map;
    std::istringstream names(functional);
    std::string name;
    while (names >> name) {
      funcs.emplace_back();
      BOOST_REQUIRE(
          xc_func_init(&funcs.back(), map.getID(name), polarization) == 0);
    }
  }
  ~Libxc() {
    for (auto& f : funcs) {
      xc_func_end(&f);
    }
  }
  // sums all functionals; arrays in libxc layout for one point
  void Eval(const double* rho, const double* sigma, double& exc, double* vrho,
            double* vsigma, int nrho, int nsigma) {
    exc = 0;
    std::fill(vrho, vrho + nrho, 0.0);
    std::fill(vsigma, vsigma + nsigma, 0.0);
    for (auto& f : funcs) {
      double e = 0;
      std::vector<double> vr(nrho, 0.0), vs(nsigma, 0.0);
      if (f.info->family == XC_FAMILY_LDA) {
        xc_lda_exc_vxc(&f, 1, rho, &e, vr.data());
      } else {
        xc_gga_exc_vxc(&f, 1, rho, sigma, &e, vr.data(), vs.data());
      }
      exc += e;
      for (int i = 0; i < nrho; ++i) {
        vrho[i] += vr[i];
      }
      for (int i = 0; i < nsigma; ++i) {
        vsigma[i] += vs[i];
      }
    }
  }
  std::vector<xc_func_type> funcs;
};

Mat_p_Energy PointByPoint(const Vxc_Grid& grid, const std::string& functional,
                          const Eigen::MatrixXd& dmat) {
  Libxc xc(functional, XC_UNPOLARIZED);
  Eigen::MatrixXd vxc = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  double energy = 0;
  for (votca::Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }
    const Eigen::MatrixXd D = 2 * box.ReadFromBigMatrix(dmat);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(D.rows(), D.cols());
    for (votca::Index p = 0; p < box.size(); ++p) {
      AOShell::AOValues ao = box.CalcAOValues(box.getGridPoints()[p]);
      Eigen::VectorXd temp = ao.values.transpose() * D;
      double rho = 0.5 * temp.dot(ao.values);
      const double w = box.getGridWeights()[p];
      if (rho * w < 1.e-20) {
        continue;
      }
      const Eigen::Vector3d g = temp.transpose() * ao.derivatives;
      double sigma = g.squaredNorm();
      double exc, vrho, vsigma;
      xc.Eval(&rho, &sigma, exc, &vrho, &vsigma, 1, 1);
      energy += w * rho * exc;
      Eigen::VectorXd a =
          w * (0.5 * vrho * ao.values + 2.0 * vsigma * (ao.derivatives * g));
      V += a * ao.values.transpose();
    }
    box.AddtoBigMatrix(vxc, V);
  }
  return Mat_p_Energy(energy, vxc + vxc.transpose());
}

std::array<Eigen::MatrixXd, 2> PointByPointSpin(const Vxc_Grid& grid,
                                                const std::string& functional,
                                                const Eigen::MatrixXd& da,
                                                const Eigen::MatrixXd& db,
                                                double& energy) {
  Libxc xc(functional, XC_POLARIZED);
  Eigen::MatrixXd va = Eigen::MatrixXd::Zero(da.rows(), da.cols());
  Eigen::MatrixXd vb = va;
  energy = 0;
  for (votca::Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }
    const Eigen::MatrixXd Da = box.ReadFromBigMatrix(da);
    const Eigen::MatrixXd Db = box.ReadFromBigMatrix(db);
    Eigen::MatrixXd Va = Eigen::MatrixXd::Zero(Da.rows(), Da.cols());
    Eigen::MatrixXd Vb = Va;
    for (votca::Index p = 0; p < box.size(); ++p) {
      AOShell::AOValues ao = box.CalcAOValues(box.getGridPoints()[p]);
      Eigen::VectorXd ta = Da * ao.values;
      Eigen::VectorXd tb = Db * ao.values;
      double rho[2] = {ao.values.dot(ta), ao.values.dot(tb)};
      const double w = box.getGridWeights()[p];
      if ((rho[0] + rho[1]) * w < 1.e-20) {
        continue;
      }
      const Eigen::Vector3d ga = 2.0 * (ao.derivatives.transpose() * ta);
      const Eigen::Vector3d gb = 2.0 * (ao.derivatives.transpose() * tb);
      double sigma[3] = {ga.dot(ga), ga.dot(gb), gb.dot(gb)};
      double exc, vrho[2], vs[3];
      xc.Eval(rho, sigma, exc, vrho, vs, 2, 3);
      energy += w * (rho[0] + rho[1]) * exc;
      Eigen::VectorXd g_a = ao.derivatives * ga;
      Eigen::VectorXd g_b = ao.derivatives * gb;
      Eigen::VectorXd wa =
          w * (0.5 * vrho[0] * ao.values + 2.0 * vs[0] * g_a + vs[1] * g_b);
      Eigen::VectorXd wb =
          w * (0.5 * vrho[1] * ao.values + vs[1] * g_a + 2.0 * vs[2] * g_b);
      Va += wa * ao.values.transpose();
      Vb += wb * ao.values.transpose();
    }
    box.AddtoBigMatrix(va, Va);
    box.AddtoBigMatrix(vb, Vb);
  }
  return {va + va.transpose(), vb + vb.transpose()};
}

void CheckClose(const Eigen::MatrixXd& result, const Eigen::MatrixXd& ref,
                const std::string& what) {
  const double err = (result - ref).cwiseAbs().maxCoeff();
  BOOST_TEST_MESSAGE(what << ": max deviation " << err);
  BOOST_CHECK_MESSAGE(err < 1e-11 * ref.cwiseAbs().maxCoeff(),
                      what << ": max deviation " << err);
}
}  // namespace

// The blocked integration must reproduce the point-by-point one to rounding,
// for a separate exchange/correlation pair, a hybrid, an LDA, and the
// spin-polarized path with unequal alpha and beta densities.
BOOST_AUTO_TEST_CASE(blocked_integration_matches_point_by_point) {
  libint2::initialize();
  QMMolecule mol("none", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/vxc_potential/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);
  const Eigen::MatrixXd dmat = DMat();
  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);

  for (const std::string functional :
       {"XC_GGA_X_PBE XC_GGA_C_PBE", "XC_HYB_GGA_XC_PBEH",
        "XC_LDA_X XC_LDA_C_VWN"}) {
    Vxc_Potential<Vxc_Grid> num(grid);
    num.setXCfunctional(functional);

    const Mat_p_Energy blocked = num.IntegrateVXC(dmat);
    const Mat_p_Energy ref = PointByPoint(grid, functional, dmat);
    CheckClose(blocked.matrix(), ref.matrix(), functional + " Vxc");
    BOOST_CHECK_SMALL(blocked.energy() - ref.energy(),
                      1e-11 * std::abs(ref.energy()));

    // unequal spin densities, and the closed-shell limit
    for (double fa : {0.6, 0.5}) {
      const Eigen::MatrixXd da = fa * dmat;
      const Eigen::MatrixXd db = (1.0 - fa) * dmat;
      const auto spin = num.IntegrateVXCSpin(da, db);
      double e_ref = 0;
      const auto ref_spin = PointByPointSpin(grid, functional, da, db, e_ref);
      CheckClose(spin.vxc_alpha, ref_spin[0], functional + " Vxc alpha");
      CheckClose(spin.vxc_beta, ref_spin[1], functional + " Vxc beta");
      BOOST_CHECK_SMALL(spin.energy - e_ref, 1e-11 * std::abs(e_ref));
      if (fa == 0.5) {
        // UKS with equal spins must equal RKS
        CheckClose(spin.vxc_alpha, blocked.matrix(),
                   functional + " UKS(equal spins) vs RKS");
        BOOST_CHECK_SMALL(spin.energy - blocked.energy(),
                          1e-10 * std::abs(blocked.energy()));
      }
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
