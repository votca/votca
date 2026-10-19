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

#define BOOST_TEST_MODULE staticresponse_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_potential.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(staticresponse_test)

namespace {

const std::string kBasis =
    std::string(XTP_TEST_DATA_FOLDER) + "/hirshfeldpartition/3-21G.xml";

QMMolecule Water() {
  QMMolecule mol("water", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/chargeresponse/water.xyz");
  return mol;
}

AOBasis Basis(const QMMolecule& mol) {
  BasisSet bs;
  bs.Load(kBasis);
  AOBasis basis;
  basis.Fill(bs, mol);
  return basis;
}

// xtpdft options; without RI if ri is false (4-centre integrals)
votca::tools::Property DFTOptions(const std::string& functional, bool ri) {
  const std::string file = "staticresponse_test.xml";
  {
    std::ofstream xml(file);
    xml << "<dftpackage>\n<spin>1</spin>\n<name>xtp</name>\n"
        << "<charge>0</charge>\n"
        << "<functional>" << functional << "</functional>\n"
        << "<basisset>" << kBasis << "</basisset>\n";
    if (ri) {
      xml << "<auxbasisset>" << XTP_TEST_DATA_FOLDER
          << "/diabatization/aux-def2-svp.xml</auxbasisset>\n";
    }
    xml << "<initial_guess>independent</initial_guess>\n"
        << "<xtpdft>\n<screening_eps>1e-12</screening_eps>\n"
        << "<fock_matrix_reset>5</fock_matrix_reset>\n<convergence>\n"
        << "<energy>1e-13</energy>\n<method>DIIS</method>\n"
        << "<DIIS_start>0.002</DIIS_start>\n<ADIIS_start>0.8</ADIIS_start>\n"
        << "<DIIS_length>20</DIIS_length>\n<levelshift>0.0</levelshift>\n"
        << "<levelshift_end>0.2</levelshift_end>\n"
        << "<max_iterations>100</max_iterations>\n<error>1e-11</error>\n"
        << "<DIIS_maxout>false</DIIS_maxout>\n<mixing>0.7</mixing>\n"
        << "<mixing_max>0.98</mixing_max>\n<mixing_end>0.8</mixing_end>\n"
        << "<davidson_max_iter>50</davidson_max_iter>\n"
        << "</convergence>\n<integration_grid>medium</integration_grid>\n"
        << "<max_iterations>200</max_iterations>\n</xtpdft>\n"
        << "</dftpackage>\n";
  }
  votca::tools::Property prop;
  prop.LoadFromXML(file);
  return prop.get("dftpackage");
}

// SCF of water with an optional electron potential energy operator v
Orbitals Ground(votca::tools::Property options, const Eigen::MatrixXd& v,
                const Orbitals* guess, Logger& log) {
  DFTEngine dft;
  dft.setLogger(&log);
  dft.Initialize(options);
  Orbitals orb;
  if (guess) {
    orb = *guess;
    dft.setWarmStart(true);
  } else {
    orb.QMAtoms() = Water();
    orb.setChargeAndSpin(0, 1);
  }
  if (v.size() > 0) {
    dft.setExternalPotentialMatrix(v);
  }
  BOOST_REQUIRE(dft.Evaluate(orb));
  return orb;
}

Eigen::MatrixXd RandomSymmetric(Index n, double scale) {
  const Eigen::MatrixXd r = Eigen::MatrixXd::Random(n, n);
  return scale * 0.5 * (r + r.transpose());
}
}  // namespace

BOOST_AUTO_TEST_CASE(fxc_is_the_derivative_of_vxc) {
  // IntegrateFXC against central differences of IntegrateVXC, for an LDA
  // and a GGA (exchange and correlation separate) and a hybrid GGA
  libint2::initialize();
  Logger log;
  const Orbitals orb = Ground(DFTOptions("XC_GGA_X_PBE XC_GGA_C_PBE", true),
                              Eigen::MatrixXd(), nullptr, log);
  const QMMolecule mol = Water();
  const AOBasis basis = Basis(mol);
  Vxc_Grid grid;
  grid.GridSetup("medium", mol, basis);
  const Eigen::MatrixXd P = orb.DensityMatrixGroundState();
  const Index n = basis.AOBasisSize();
  // a first-order density with zero trace against S: no change of the
  // electron number, like the CPKS densities
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::MatrixXd dP = RandomSymmetric(n, 0.05);
  dP -= (dP.cwiseProduct(overlap.Matrix()).sum() /
         P.cwiseProduct(overlap.Matrix()).sum()) *
        P;
  for (const std::string functional :
       {"XC_LDA_X XC_LDA_C_VWN", "XC_GGA_X_PBE XC_GGA_C_PBE",
        "XC_HYB_GGA_XC_PBEH"}) {
    Vxc_Potential<Vxc_Grid> vxc(grid);
    vxc.setXCfunctional(functional);
    const Eigen::MatrixXd fxc = vxc.IntegrateFXC(P, {dP})[0];
    const double h = 1e-4;
    const Eigen::MatrixXd fd = (vxc.IntegrateVXC(P + h * dP).matrix() -
                                vxc.IntegrateVXC(P - h * dP).matrix()) /
                               (2 * h);
    const double err = (fxc - fd).cwiseAbs().maxCoeff();
    BOOST_TEST_MESSAGE(functional << ": max |fxc - fd| " << err << " of "
                                  << fd.cwiseAbs().maxCoeff());
    BOOST_CHECK_SMALL(err, 1e-6 * fd.cwiseAbs().maxCoeff());
    BOOST_CHECK_SMALL((fxc - fxc.transpose()).cwiseAbs().maxCoeff(), 1e-14);
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(response_matches_perturbed_scfs) {
  // first-order densities from CPKS against central differences of SCFs in
  // the same perturbations (a uniform field and a localized potential),
  // for a hybrid with RI and a GGA with 4-centre integrals
  libint2::initialize();
  Logger log;
  const QMMolecule mol = Water();
  const AOBasis basis = Basis(mol);
  AODipole dipole;
  dipole.setCenter(mol.getPos());
  dipole.Fill(basis);
  // a Gaussian-weighted potential on the oxygen: diagonal-dominant, smooth
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::MatrixXd local =
      Eigen::MatrixXd::Zero(basis.AOBasisSize(), basis.AOBasisSize());
  for (const AOShell& shell : basis) {
    if (shell.getAtomIndex() == 0) {
      const Index s = shell.getStartIndex();
      for (Index f = 0; f < shell.getNumFunc(); ++f) {
        local(s + f, s + f) = 1.0;
      }
    }
  }
  local = 0.5 * (local * overlap.Matrix() + overlap.Matrix() * local);
  const std::vector<Eigen::MatrixXd> h = {dipole.Matrix()[2], local};

  for (const bool hybrid_ri : {true, false}) {
    const votca::tools::Property options =
        hybrid_ri ? DFTOptions("XC_HYB_GGA_XC_PBEH", true)
                  : DFTOptions("XC_GGA_X_PBE XC_GGA_C_PBE", false);
    const Orbitals orb = Ground(options, Eigen::MatrixXd(), nullptr, log);
    DFTEngine dft;
    dft.setLogger(&log);
    votca::tools::Property opt = options;
    dft.Initialize(opt);
    const std::vector<Eigen::MatrixXd> dP = dft.StaticResponse(orb, h, 1e-9);
    BOOST_REQUIRE_EQUAL(dP.size(), 2);
    const double step = 1e-4;
    for (std::size_t k = 0; k < h.size(); ++k) {
      const Orbitals plus = Ground(options, step * h[k], &orb, log);
      const Orbitals minus = Ground(options, -step * h[k], &orb, log);
      const Eigen::MatrixXd fd =
          (plus.DensityMatrixGroundState() - minus.DensityMatrixGroundState()) /
          (2 * step);
      const double err = (dP[k] - fd).cwiseAbs().maxCoeff();
      BOOST_TEST_MESSAGE((hybrid_ri ? "PBE0/RI" : "PBE/4c")
                         << " perturbation " << k << ": max |dP - fd| " << err
                         << " of " << fd.cwiseAbs().maxCoeff());
      BOOST_CHECK_SMALL(err, 1e-5 * fd.cwiseAbs().maxCoeff());
      // the electron number does not change
      BOOST_CHECK_SMALL(dP[k].cwiseProduct(overlap.Matrix()).sum(), 1e-10);
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
