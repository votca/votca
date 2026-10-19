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
#include <filesystem>
#include <fstream>
#include <iomanip>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/orbitals.h"

// Local private VOTCA includes
#include "crfit.h"

namespace votca {
namespace xtp {

void CRFit::ParseOptions(const tools::Property& options) {
  xyzfile_ = job_name_ + ".xyz";
  dftoptions_ = options.get(".dftpackage");
  if (dftoptions_.get("name").as<std::string>() != "xtp") {
    throw std::runtime_error("crfit: needs the xtp dftpackage");
  }
  step_ = options.get(".step").as<double>();
  scf_error_ = options.get(".scf_error").as<double>();
  inversion_ = options.get(".inversion").as<bool>();
  output_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".output", job_name_ + "_cr");
  if (options.exists(".mpsfile")) {
    mpsfile_ = options.get(".mpsfile").as<std::string>();
  }
  // each SCF at least as tight as scf_error (the differences divide by step)
  tools::Property& conv = dftoptions_.get("xtpdft.convergence");
  const double error = conv.get("error").as<double>();
  auto str = [](double x) {
    std::stringstream s;
    s << std::setprecision(17) << x;
    return s.str();
  };
  if (error > scf_error_) {
    conv.set("error", str(scf_error_));
  }
  const double energy = conv.get("energy").as<double>();
  if (energy > 1e-2 * scf_error_) {
    conv.set("energy", str(1e-2 * scf_error_));
  }
}

bool CRFit::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  Orbitals orb0;
  orb0.QMAtoms().LoadFromFile(xyzfile_);
  const Index charge = dftoptions_.get("charge").as<Index>();
  const Index spin = dftoptions_.get("spin").as<Index>();
  if (spin != 1) {
    throw std::runtime_error(
        "crfit: closed-shell molecules only (spin 1), got spin " +
        std::to_string(spin));
  }
  orb0.setChargeAndSpin(charge, spin);
  const QMMolecule& mol = orb0.QMAtoms();
  const Index n_atoms = mol.size();

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Charge-response fit of " << n_atoms << " atoms from "
      << xyzfile_ << std::flush;
  const ChargeResponseFit::Result r =
      ChargeResponseFit::Fit(dftoptions_, orb0, step_, log_, inversion_);

  const double A3 = std::pow(tools::conv::bohr2ang, 3);
  auto eig = [&](const Eigen::Matrix3d& m) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(m,
                                                      Eigen::EigenvaluesOnly);
    return Eigen::Vector3d(es.eigenvalues() * A3);
  };
  XTP_LOG(Log::error, log_)
      << std::setprecision(4) << TimeStamp()
      << " alpha eigenvalues [A^3]: " << eig(r.alpha).transpose()
      << "\n charge flow X^T Pi X: " << eig(r.alpha_ct).transpose()
      << "\n residual: " << eig(r.alpha_res).transpose()
      << "\n Pi: asymmetry before symmetrizing " << r.pi_asymmetry
      << ", eigenvalues " << r.pi_eigenvalues(1) << " .. "
      << r.pi_eigenvalues(n_atoms - 1) << " e^2/Hartree" << std::flush;

  WriteCrm(mol, r);
  WriteJson(mol, r);
  if (!mpsfile_.empty()) {
    WriteMps(mol, r);
  }
  return true;
}

void CRFit::WriteCrm(const QMMolecule& mol,
                     const ChargeResponseFit::Result& r) const {
  const std::string file = output_ + ".crm";
  std::ofstream ofs(file);
  if (!ofs.is_open()) {
    throw std::runtime_error("crfit: cannot write " + file);
  }
  ofs << "! charge-response matrix (Hirshfeld), "
      << dftoptions_.get("functional").as<std::string>() << "/"
      << dftoptions_.get("basisset").as<std::string>() << ", from " << xyzfile_
      << " (xtp_tools -e crfit, finite differences, step " << step_ << ")\n";
  ofs << "! dq = -Pi phi; Pi in e^2/Hartree, atom order as below (Angstrom)\n";
  ofs << mol.size() << "\n";
  for (const QMAtom& atom : mol) {
    const Eigen::Vector3d p = atom.getPos() * tools::conv::bohr2ang;
    ofs << boost::format("%-2s %12.6f %12.6f %12.6f\n") % atom.getElement() %
               p(0) % p(1) % p(2);
  }
  for (Index j = 0; j < r.pi.rows(); ++j) {
    for (Index k = 0; k < r.pi.cols(); ++k) {
      ofs << boost::format("%+.10e") % r.pi(j, k)
          << (k + 1 < r.pi.cols() ? " " : "\n");
    }
  }
  XTP_LOG(Log::error, log_) << TimeStamp() << " Wrote " << file << std::flush;
}

namespace {
std::string Json(const Eigen::MatrixXd& m) {
  std::stringstream s;
  s << std::setprecision(15) << "[";
  for (Index i = 0; i < m.rows(); ++i) {
    s << (i ? ", [" : "[");
    for (Index j = 0; j < m.cols(); ++j) {
      s << (j ? ", " : "") << m(i, j);
    }
    s << "]";
  }
  return s.str() + "]";
}
std::string JsonVector(const Eigen::VectorXd& v) {
  std::stringstream s;
  s << std::setprecision(15);
  for (Index i = 0; i < v.size(); ++i) {
    s << (i ? ", " : "") << v(i);
  }
  return s.str();
}
}  // namespace

void CRFit::WriteJson(const QMMolecule& mol,
                      const ChargeResponseFit::Result& r) const {
  // the keys of xtp_cr_fit, so that xtp_cr_mps reads both
  const std::string file = output_ + ".json";
  std::ofstream ofs(file);
  if (!ofs.is_open()) {
    throw std::runtime_error("crfit: cannot write " + file);
  }
  const double A3 = std::pow(tools::conv::bohr2ang, 3);
  Eigen::MatrixXd coords(mol.size(), 3);
  std::string elements;
  for (Index k = 0; k < mol.size(); ++k) {
    coords.row(k) = mol[k].getPos().transpose();
    elements += (k ? ", \"" : "\"") + mol[k].getElement() + "\"";
  }
  std::string res_atoms;
  for (std::size_t k = 0; k < r.alpha_res_atoms.size(); ++k) {
    res_atoms += (k ? ", " : "") + Json(r.alpha_res_atoms[k]);
  }
  ofs << "{\n"
      << " \"xyz\": \"" << xyzfile_ << "\",\n"
      << " \"xc\": \"" << dftoptions_.get("functional").as<std::string>()
      << "\",\n"
      << " \"basis\": \"" << dftoptions_.get("basisset").as<std::string>()
      << "\",\n"
      << " \"charge\": " << dftoptions_.get("charge").as<Index>() << ",\n"
      << " \"method\": \"xtp_tools -e crfit (finite differences)\",\n"
      << " \"step\": " << step_ << ",\n"
      << " \"units\": \"atomic (bohr, e, Hartree); Pi in e^2/Hartree\",\n"
      << " \"elements\": [" << elements << "],\n"
      << " \"coords_bohr\": " << Json(coords) << ",\n"
      << " \"alpha_dft_au\": " << Json(r.alpha) << ",\n"
      << " \"alpha_dft_A3\": " << Json(r.alpha * A3) << ",\n"
      << " \"Pi\": " << Json(r.pi) << ",\n"
      << " \"Pi_eigenvalues\": [" << JsonVector(r.pi_eigenvalues) << "],\n"
      << std::setprecision(6) << " \"Pi_asymmetry\": " << r.pi_asymmetry
      << ",\n"
      << " \"alpha_ct_A3\": " << Json(r.alpha_ct * A3) << ",\n"
      << " \"alpha_res_A3\": " << Json(r.alpha_res * A3) << ",\n"
      << " \"alpha_res_weights\": [" << JsonVector(r.res_weights) << "],\n"
      << " \"alpha_res_atoms_au\": [" << res_atoms << "]\n"
      << "}\n";
  XTP_LOG(Log::error, log_) << TimeStamp() << " Wrote " << file << std::flush;
}

void CRFit::WriteMps(const QMMolecule& mol,
                     const ChargeResponseFit::Result& r) const {
  PolarSegment tmpl("template", 0);
  tmpl.LoadFromFile(mpsfile_);
  if (tmpl.size() != mol.size()) {
    throw std::runtime_error(
        "crfit: " + mpsfile_ + " has " + std::to_string(tmpl.size()) +
        " sites, the molecule " + std::to_string(mol.size()) + " atoms");
  }
  Eigen::MatrixX3d P(mol.size(), 3), Q(mol.size(), 3);
  for (Index k = 0; k < mol.size(); ++k) {
    if (tmpl[k].getElement() != mol[k].getElement()) {
      throw std::runtime_error("crfit: site " + std::to_string(k) + " of " +
                               mpsfile_ + " is " + tmpl[k].getElement() +
                               ", atom " + std::to_string(k) + " is " +
                               mol[k].getElement());
    }
    P.row(k) = tmpl[k].getPos().transpose();
    Q.row(k) = mol[k].getPos().transpose();
  }
  // rigid alignment (Kabsch) of the fit geometry onto the template, to
  // rotate the residual tensors into the template frame
  const Eigen::MatrixX3d P0 = P.rowwise() - P.colwise().mean();
  const Eigen::MatrixX3d Q0 = Q.rowwise() - Q.colwise().mean();
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(
      Q0.transpose() * P0, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
  D(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant() < 0 ? -1.0
                                                                          : 1.0;
  const Eigen::Matrix3d U =
      svd.matrixV() * D * svd.matrixU().transpose();  // fit -> template
  const double rmsd =
      std::sqrt((Q0 * U.transpose() - P0).rowwise().squaredNorm().mean());
  if (rmsd * tools::conv::bohr2ang > 0.05) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " WARNING: template and fit geometry differ by "
        << rmsd * tools::conv::bohr2ang << " A rms after alignment"
        << std::flush;
  }
  for (Index k = 0; k < mol.size(); ++k) {
    tmpl[k].setpolarization(U * r.alpha_res_atoms[std::size_t(k)] *
                            U.transpose());
  }
  const std::string file = output_ + ".mps";
  tmpl.WriteMPS(file, "CRFIT (permanent multipoles of " +
                          std::filesystem::path(mpsfile_).filename().string() +
                          ", residual polarizabilities of the fit)");
  XTP_LOG(Log::error, log_) << TimeStamp() << " Wrote " << file << std::flush;
}

}  // namespace xtp
}  // namespace votca
