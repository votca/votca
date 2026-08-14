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

  // Same index-list syntax as CDFT's own cdft.indices and
  // diabatization.xml's own fragment indices (space-separated
  // indices and/or ranges, e.g. "1 3 13:17") -- deliberately reusing
  // the SAME, already-established IndexParser utility both of those
  // already use, rather than introducing yet another parsing
  // convention for what is, at its core, the exact same kind of
  // input (an atom-index list defining a fragment).
  std::string fragment_A_str = options.get(".fragment_A").as<std::string>();
  std::string fragment_B_str = options.get(".fragment_B").as<std::string>();
  if (fragment_A_str.empty() || fragment_B_str.empty()) {
    throw std::runtime_error(
        "PodCouplingTool: fragment_A and fragment_B must both be given, "
        "e.g. '0:5' and '6:11' -- same syntax as cdft.indices and "
        "diabatization.xml's own fragment indices.");
  }
  fragment_A_atoms_ = IndexParser().CreateIndexVector(fragment_A_str);
  fragment_B_atoms_ = IndexParser().CreateIndexVector(fragment_B_str);

  // Matching DFTcoupling's own numberofstatesA_/numberofstatesB_
  // convention exactly (same XML option names, levA/levB) -- covers
  // BOTH occupied and virtual orbitals for each fragment in a single
  // call, per direct agreement with the user to mirror DFTcoupling's
  // own, established behavior here rather than this tool's own,
  // earlier, single-orbital-pair-only interface.
  numberofstatesA_ =
      options.ifExistsReturnElseReturnDefault<Index>("levA", 1);
  numberofstatesB_ =
      options.ifExistsReturnElseReturnDefault<Index>("levB", 1);

  XTP_LOG(Log::error, log_)
      << "Fragment A: " << fragment_A_atoms_.size()
      << " atoms, levA=" << numberofstatesA_ << flush;
  XTP_LOG(Log::error, log_)
      << "Fragment B: " << fragment_B_atoms_.size()
      << " atoms, levB=" << numberofstatesB_ << flush;
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

  // The full, pairwise coupling matrix across the requested range,
  // covering both hole (occupied) and electron (virtual) orbitals for
  // each fragment -- matching DFTcoupling's own output style of
  // reporting the whole JAB matrix, not just a single pair, per direct
  // agreement with the user. Rows/columns are labeled by each
  // orbital's own, absolute (fragment-local) index -- see
  // PODCoupling::getFragmentAHomoIndex/etc.'s own header comment for
  // what "absolute" means here.
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Full coupling matrix [eV] (fragment A orbital "
                        "index -> fragment B orbital index):"
      << flush;
  for (Index levelA = homoA - numberofstatesA_ + 1;
      levelA <= lumoA + numberofstatesA_ - 1; ++levelA) {
    std::ostringstream row;
    row << "  A[" << levelA << (levelA == homoA ? "=HOMO" : "")
        << (levelA == lumoA ? "=LUMO" : "") << "]:";
    for (Index levelB = homoB - numberofstatesB_ + 1;
        levelB <= lumoB + numberofstatesB_ - 1; ++levelB) {
      double coupling_ev = pod.getCouplingElement(levelA, levelB) *
                           votca::tools::conv::hrt2ev;
      row << boost::format("  B[%1%]=%2$+1.6f") % levelB % coupling_ev;
    }
    XTP_LOG(Log::error, log_) << row.str() << flush;
  }

  // The two, most commonly-needed values highlighted directly --
  // matching the paper's own default focus on hole transport
  // (HOMO-HOMO), plus its own, direct electron-transport analog
  // (LUMO-LUMO) -- so a caller does not need to parse the full matrix
  // above just to get the single, most likely value of interest.
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

  std::cout << log_;
  // Explicit trailing newline: confirmed directly, from the user's
  // own real run, that xtp_tools.cc's own execute() ends right after
  // tool_->Evaluate() with no final newline anywhere in the shared
  // framework -- the shell's own next prompt was running directly
  // into this tool's own last line of output as a result. This is a
  // local, low-risk fix within this tool's own Run() specifically,
  // NOT a fix to the shared framework itself (which affects every
  // existing tool, not just this one, and would need its own,
  // separate, deliberate change and testing before touching).
  std::cout << std::endl;

  return true;
}

}  // namespace xtp
}  // namespace votca
