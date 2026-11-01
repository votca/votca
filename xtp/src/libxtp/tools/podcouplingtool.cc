/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#include "podcouplingtool.h"
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/cubefile_writer.h"
#include "votca/xtp/qmstate.h"
#include <sstream>

using std::flush;

namespace votca {
namespace xtp {

void PodCouplingTool::ParseOptions(const tools::Property& user_options) {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n...");

  tools::Property options = user_options;

  orb_file_ = options.get(".orb_file").as<std::string>();

  // Atom index lists in IndexParser syntax (e.g. "1 3 13:17"), as for
  // cdft.indices and the fragments in diabatization.xml.
  std::string fragment_A_str = options.get(".fragment_A").as<std::string>();
  std::string fragment_B_str = options.get(".fragment_B").as<std::string>();
  if (fragment_A_str.empty() || fragment_B_str.empty()) {
    throw std::runtime_error(
        "PodCouplingTool: fragment_A and fragment_B must both be given, "
        "e.g. '0:5' and '6:11' (same syntax as cdft.indices).");
  }
  fragment_A_atoms_ = IndexParser().CreateIndexVector(fragment_A_str);
  fragment_B_atoms_ = IndexParser().CreateIndexVector(fragment_B_str);

  // Number of occupied and virtual levels per fragment; option names and
  // meaning as in DFTcoupling (levA/levB).
  numberofstatesA_ = options.ifExistsReturnElseReturnDefault<Index>("levA", 1);
  numberofstatesB_ = options.ifExistsReturnElseReturnDefault<Index>("levB", 1);

  XTP_LOG(Log::error, log_) << "Fragment A: " << fragment_A_atoms_.size()
                            << " atoms, levA=" << numberofstatesA_ << flush;
  XTP_LOG(Log::error, log_) << "Fragment B: " << fragment_B_atoms_.size()
                            << " atoms, levB=" << numberofstatesB_ << flush;

  // Optional cube-file export; grid options as in gencube.xml.
  write_cube_files_ =
      options.ifExistsReturnElseReturnDefault<bool>("write_cube_files", false);
  if (write_cube_files_) {
    cube_steps_.x() =
        options.ifExistsReturnElseReturnDefault<Index>("xsteps", 25);
    cube_steps_.y() =
        options.ifExistsReturnElseReturnDefault<Index>("ysteps", 25);
    cube_steps_.z() =
        options.ifExistsReturnElseReturnDefault<Index>("zsteps", 25);
    cube_padding_ =
        options.ifExistsReturnElseReturnDefault<double>("padding", 6.5);
  }
}

bool PodCouplingTool::Run() {

  OPENMP::setMaxThreads(nThreads_);

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n...");

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Reading orbitals from: " << orb_file_ << flush;

  Orbitals orbitals;
  orbitals.ReadFromCpt(orb_file_);

  PODCoupling pod(orbitals, &log_, fragment_A_atoms_, fragment_B_atoms_);
  pod.CalculateCouplings(numberofstatesA_, numberofstatesB_);

  Index homoA = pod.getFragmentAHomoIndex();
  Index lumoA = pod.getFragmentALumoIndex();
  Index homoB = pod.getFragmentBHomoIndex();
  Index lumoB = pod.getFragmentBLumoIndex();

  // Full coupling matrix over the requested range, as in DFTcoupling,
  // labeled by absolute fragment orbital indices.
  XTP_LOG(Log::error, log_) << TimeStamp()
                            << " Full coupling matrix [eV] (fragment A orbital "
                               "index -> fragment B orbital index):"
                            << flush;
  for (Index levelA = homoA - numberofstatesA_ + 1;
       levelA <= lumoA + numberofstatesA_ - 1; ++levelA) {
    std::ostringstream row;
    row << "  A[" << levelA << (levelA == homoA ? "=HOMO" : "")
        << (levelA == lumoA ? "=LUMO" : "") << "]:";
    for (Index levelB = homoB - numberofstatesB_ + 1;
         levelB <= lumoB + numberofstatesB_ - 1; ++levelB) {
      double coupling_ev =
          pod.getCouplingElement(levelA, levelB) * votca::tools::conv::hrt2ev;
      row << boost::format("  B[%1%]=%2$+1.6f") % levelB % coupling_ev;
    }
    XTP_LOG(Log::error, log_) << row.str() << flush;
  }

  // HOMO-HOMO (hole) and LUMO-LUMO (electron) couplings separately.
  double homo_homo_ev =
      pod.getCouplingElement(homoA, homoB) * votca::tools::conv::hrt2ev;
  double lumo_lumo_ev =
      pod.getCouplingElement(lumoA, lumoB) * votca::tools::conv::hrt2ev;
  XTP_LOG(Log::error, log_)
      << boost::format("POD2 Coupling, HOMO-HOMO (hole): %1$+1.12f eV") %
             homo_homo_ev
      << flush;
  XTP_LOG(Log::error, log_)
      << boost::format("POD2 Coupling, LUMO-LUMO (electron): %1$+1.12f eV") %
             lumo_lumo_ev
      << flush;

  // Cube files of the computed fragment orbitals. CubeFile_Writer only
  // reads orb.MOs().eigenvectors().col(index), so each orbital, embedded
  // in the full AO basis (PODCoupling::GetFragmentOrbital), is stored as
  // the single MO of a copy of orbitals and written as in gencube.cc.
  if (write_cube_files_) {
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " Writing cube files for the computed fragment "
           "orbitals..."
        << flush;
    CubeFile_Writer writer(cube_steps_, cube_padding_, log_);
    auto WriteFragmentOrbitalCube = [&](bool fragment_A, Index level,
                                        const std::string& label) {
      Eigen::VectorXd orbital = pod.GetFragmentOrbital(fragment_A, level);
      Orbitals orbitals_copy = orbitals;
      orbitals_copy.MOs().eigenvectors() = orbital;
      orbitals_copy.MOs().eigenvalues() = Eigen::VectorXd::Zero(1);
      std::string filename =
          "podcoupling_" + label + "_orb" + std::to_string(level) + ".cube";
      writer.WriteFile(filename, orbitals_copy,
                       QMState(QMStateType::KSstate, 0, false), false);
      XTP_LOG(Log::error, log_) << "  wrote " << filename << flush;
      // Dominant AOs of the orbital, e.g. to identify which AO's tail
      // causes density in the region of the other fragment.
      XTP_LOG(Log::error, log_)
          << pod.DescribeFragmentOrbitalComposition(fragment_A, level) << flush;
    };
    for (Index levelA = homoA - numberofstatesA_ + 1;
         levelA <= lumoA + numberofstatesA_ - 1; ++levelA) {
      WriteFragmentOrbitalCube(true, levelA, "fragA");
    }
    for (Index levelB = homoB - numberofstatesB_ + 1;
         levelB <= lumoB + numberofstatesB_ - 1; ++levelB) {
      WriteFragmentOrbitalCube(false, levelB, "fragB");
    }
  }

  std::cout << log_;
  // Trailing newline: xtp_tools does not print one after Evaluate(), so
  // the shell prompt would follow the last output line directly.
  std::cout << std::endl;

  return true;
}

}  // namespace xtp
}  // namespace votca
