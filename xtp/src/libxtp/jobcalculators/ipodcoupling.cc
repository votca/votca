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
#include <limits>
#include <set>

// Third party includes
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/xtp/fragmentsaturator.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/podcoupling.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/segmentmapper.h"

// Local private VOTCA includes
#include "ipodcoupling.h"

namespace votca {
namespace xtp {

void IPodCoupling::ParseSpecificOptions(const tools::Property& options) {
  QMPackageFactory{};

  // Task names follow IQM::ParseSpecificOptions (iqm.cc), plus the
  // additional "podcoupling" task.
  std::string tasks_string = options.get(".tasks").as<std::string>();

  // We split either on a space or a comma
  tools::Tokenizer tokenizedTasks(tasks_string, " ,");
  std::vector<std::string> tasks = tokenizedTasks.ToVector();

  do_dft_input_ = std::find(tasks.begin(), tasks.end(), "input") != tasks.end();
  do_dft_run_ = std::find(tasks.begin(), tasks.end(), "dft") != tasks.end();
  do_dft_parse_ = std::find(tasks.begin(), tasks.end(), "parse") != tasks.end();
  do_podcoupling_ =
      std::find(tasks.begin(), tasks.end(), "podcoupling") != tasks.end();

  store_dft_ =
      options.ifExistsReturnElseReturnDefault<bool>(".store_dft", store_dft_);
  include_linkers_ = options.ifExistsReturnElseReturnDefault<bool>(
      ".include_linkers", include_linkers_);

  dftpackage_options_ = options.get(".dftpackage");
  if (options.exists(".podcoupling")) {
    podcoupling_options_ = options.get(".podcoupling");
    numberofstatesA_ =
        podcoupling_options_.ifExistsReturnElseReturnDefault<Index>(
            "levA", numberofstatesA_);
    numberofstatesB_ =
        podcoupling_options_.ifExistsReturnElseReturnDefault<Index>(
            "levB", numberofstatesB_);
  }
}

void IPodCoupling::WriteJobFile(const Topology& top) {
  // One job per neighbor-list pair, as in IQM::WriteJobFile (iqm.cc).
  std::cout << std::endl
            << "... ... Writing job file " << jobfile_ << std::flush;
  std::ofstream ofs;
  ofs.open(jobfile_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("\nERROR: bad file handle: " + jobfile_);
  }

  const QMNBList& nblist = top.NBList();

  Index jobCount = 0;
  if (nblist.size() == 0) {
    std::cout << std::endl
              << "... ... No pairs in neighbor list, skip." << std::flush;
    return;
  }

  ofs << "<jobs>" << std::endl;
  std::string tag = "";

  for (const QMPair* pair : nblist) {
    if (pair->getType() == QMPair::Excitoncl) {
      continue;
    }
    Index id1 = pair->Seg1()->getId();
    std::string name1 = pair->Seg1()->getType();
    Index id2 = pair->Seg2()->getId();
    std::string name2 = pair->Seg2()->getType();
    Index id = jobCount;
    tools::Property Input;
    tools::Property& pInput = Input.add("input", "");
    tools::Property& pSegmentA =
        pInput.add("segment", boost::lexical_cast<std::string>(id1));
    pSegmentA.setAttribute<std::string>("type", name1);
    pSegmentA.setAttribute<Index>("id", id1);
    tools::Property& pSegmentB =
        pInput.add("segment", boost::lexical_cast<std::string>(id2));
    pSegmentB.setAttribute<std::string>("type", name2);
    pSegmentB.setAttribute<Index>("id", id2);
    Job job(id, tag, Input, Job::AVAILABLE);
    job.ToStream(ofs);
    jobCount++;
  }
  ofs << "</jobs>" << std::endl;
  ofs.close();
  std::cout << std::endl
            << "... ... In total " << jobCount << " jobs" << std::flush;
  return;
}

const Atom& IPodCoupling::FindBoundaryAtomTowardSegment(
    const Segment& seg, Index target_segment_id) const {
  for (const Atom& atom : seg) {
    if (atom.hasExternalBond() &&
        atom.getExternalBondPartnerSegmentId() == target_segment_id) {
      return atom;
    }
  }
  throw std::runtime_error(
      "IPodCoupling::FindBoundaryAtomTowardSegment: segment " +
      std::to_string(seg.getId()) +
      " has no atom with an external bond toward segment " +
      std::to_string(target_segment_id) +
      ", although FindLinkingSegments() reported a bond path between them. "
      "Check the underlying mapping/checkpoint data.");
}

std::vector<Segment> IPodCoupling::PositionLinkersAlongChain(
    const Topology& top, const Segment& seg1_positioned,
    const std::vector<const Segment*>& linkers,
    const Segment& seg2_positioned) const {
  std::vector<Segment> results;
  if (linkers.empty()) {
    return results;
  }
  // Required, not an optimization: previous_positioned points into
  // results, and a reallocating push_back() would invalidate it.
  results.reserve(linkers.size());

  // Walk the chain one covalent bond at a time instead of applying a
  // single whole-chain shift: each linker may be wrapped into a different
  // periodic image than its neighbors. Each hop is a bond length, far
  // shorter than the box, so the minimum-image convention is unambiguous
  // at every step.
  const Segment* previous_positioned = &seg1_positioned;
  for (const Segment* current_raw : linkers) {
    const Atom& prev_boundary = FindBoundaryAtomTowardSegment(
        *previous_positioned, current_raw->getId());
    const Atom& curr_boundary = FindBoundaryAtomTowardSegment(
        *current_raw, previous_positioned->getId());

    // Approximate position of the current segment's boundary atom, seen
    // from the already positioned previous segment's side of the bond
    // (same construction as FragmentSaturator::SaturateExternalBonds).
    Eigen::Vector3d expected_partner_pos =
        prev_boundary.getPos() +
        FragmentSaturator::kDefaultCHBondLengthAngstrom *
            tools::conv::ang2bohr * prev_boundary.getExternalBondDirection();

    // Shift for the whole current segment that moves its boundary atom
    // to the closest periodic image of expected_partner_pos.
    // PbShortestConnect(r_i, r_j) returns r_j - r_i (minimum image), so
    // the argument order matters.
    Eigen::Vector3d shift =
        top.PbShortestConnect(curr_boundary.getPos(), expected_partner_pos);

    Segment current_shifted = *current_raw;
    current_shifted.Translate(shift);
    results.push_back(current_shifted);
    previous_positioned = &results.back();
  }

  // Validate the last hop into seg2_positioned as well (throws if the
  // last linker is not bonded to it). seg2_positioned itself is not
  // shifted; it is already positioned, typically via QMPair::Seg2PbCopy().
  FindBoundaryAtomTowardSegment(*previous_positioned, seg2_positioned.getId());
  FindBoundaryAtomTowardSegment(seg2_positioned, previous_positioned->getId());

  return results;
}

Job::JobResult IPodCoupling::EvalJob(const Topology& top, Job& job,
                                     QMThread& opThread) {
  Job::JobResult jres = Job::JobResult();
  Logger& pLog = opThread.getLogger();

  std::string ipodcoupling_work_dir = "OR_FILES";
  std::string frame_dir =
      "frame_" + boost::lexical_cast<std::string>(top.getStep());

  QMMapper mapper(pLog);
  mapper.LoadMappingFile(mapfile_);

  // Segment ids and optional "qm_geometry" states (default: ground
  // state "n"), as in IQM::EvalJob (iqm.cc).
  tools::Property job_input = job.getInput();
  std::vector<tools::Property*> segment_list = job_input.Select("segment");
  Index ID_A = segment_list.front()->getAttribute<Index>("id");
  Index ID_B = segment_list.back()->getAttribute<Index>("id");

  std::string qmgeo_state_A = "n";
  if (segment_list.front()->exists("qm_geometry")) {
    qmgeo_state_A =
        segment_list.front()->getAttribute<std::string>("qm_geometry");
  }
  std::string qmgeo_state_B = "n";
  if (segment_list.back()->exists("qm_geometry")) {
    qmgeo_state_B =
        segment_list.back()->getAttribute<std::string>("qm_geometry");
  }
  QMState stateA(qmgeo_state_A);
  QMState stateB(qmgeo_state_B);

  const Segment& seg_A = top.getSegment(ID_A);
  const Segment& seg_B = top.getSegment(ID_B);
  const QMNBList& nblist = top.NBList();
  const QMPair* pair = nblist.FindPair(&seg_A, &seg_B);
  if (pair == nullptr) {
    SetJobToFailed(jres, pLog,
                   "No pair " + std::to_string(ID_A) + ":" +
                       std::to_string(ID_B) + " found in the neighbor list.");
    return jres;
  }

  // Path naming as in IQM::EvalJob (iqm.cc), but without the monomer
  // orbital files (orbFileA/orbFileB, eqm_work_dir): there is no dimer
  // guess, since monomer calculations are not wanted or not possible here.
  // "pairs_ipodcoupling" avoids collisions with IQM's "pairs_iqm".
  std::string pair_dir =
      (boost::format("%1%%2%%3%%4%%5%") % "pair" % "_" % ID_A % "_" % ID_B)
          .str();
  std::filesystem::path arg_path;
  std::string orbFileAB =
      (arg_path / ipodcoupling_work_dir / "pairs_ipodcoupling" / frame_dir /
       (boost::format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B % ".orb")
           .str())
          .generic_string();
  std::string package_append = "workdir_" + Identify();
  std::string work_dir =
      (arg_path / ipodcoupling_work_dir / package_append / frame_dir / pair_dir)
          .generic_string();

  // PBC-correct positioning of the pair via QMPair::Seg2PbCopy(), as in
  // IQM::EvalJob.
  const Segment* seg1 = pair->Seg1();
  Segment seg2 = pair->Seg2PbCopy();

  // Repeats the check of Md2QmEngine::map (md2qmengine.cc), whose warning
  // is printed only once at xtp_map time and is easily missed: without
  // MD-level bond data (Atom::getBondedPartnerIds()) no external bonds are
  // detected, so H-saturation does nothing and dangling valences reach the
  // DFT calculation.
  bool seg1_has_any_bonds = false;
  for (const Atom& atom : *seg1) {
    if (atom.getBondedPartnerIds()[0] != -1) {
      seg1_has_any_bonds = true;
      break;
    }
  }
  bool seg2_has_any_bonds = false;
  for (const Atom& atom : seg2) {
    if (atom.getBondedPartnerIds()[0] != -1) {
      seg2_has_any_bonds = true;
      break;
    }
  }
  if (!seg1_has_any_bonds && !seg2_has_any_bonds) {
    XTP_LOG(Log::error, pLog)
        << "WARNING: neither segment " << seg1->getId() << " nor segment "
        << seg_B.getId()
        << " has any MD-level bond connectivity; H-saturation of cut "
           "segment boundaries will do nothing for this pair. The topology "
           "reader used at xtp_map time probably provided no bond data."
        << std::flush;
  }

  // Linker segments: discovery (Topology::FindLinkingSegments) and
  // PBC-correct positioning (PositionLinkersAlongChain). Declared before
  // segments, which holds pointers into positioned_linkers.
  std::vector<Segment> positioned_linkers;
  if (include_linkers_) {
    std::vector<const Segment*> linkers = top.FindLinkingSegments(*seg1, seg_B);
    if (!linkers.empty()) {
      positioned_linkers = PositionLinkersAlongChain(top, *seg1, linkers, seg2);
    }
  }

  std::vector<const Segment*> segments = {seg1, &seg2};
  for (const Segment& linker : positioned_linkers) {
    segments.push_back(&linker);
  }

  QMMolecule qmmol = mapper.map(*seg1, stateA);
  // Mapped atom count of fragment A. Not taken from seg1->size(): the
  // mapping need not map every MD atom.
  Index n_fragment_A_atoms = qmmol.size();
  qmmol.AddContainer(mapper.map(seg2, stateB));
  // Mapped atom count of fragments A and B without linker atoms (unlike
  // n_original_atoms below); upper bound of fragment B's atom range.
  Index n_fragment_AB_atoms = qmmol.size();

  // Linkers are always mapped in the ground state ("n"): they are
  // selected by bond connectivity, not by segment type, so no per-type
  // state map as in IQM's linker_names is needed, and they are treated
  // as neutral bridges not involved in the charge transfer.
  //
  // linker_atom_ids holds the pre-saturation atom ids of all linker
  // atoms; these belong to neither fragment_A_atoms nor fragment_B_atoms.
  std::set<Index> linker_atom_ids;
  for (const Segment& linker : positioned_linkers) {
    Index linker_start_atom_id = qmmol.size();
    qmmol.AddContainer(mapper.map(linker, QMState("n")));
    for (Index i = linker_start_atom_id; i < qmmol.size(); i++) {
      linker_atom_ids.insert(i);
    }
  }

  // An external bond whose partner segment is part of this supermolecule
  // is satisfied internally and must not be saturated with H.
  std::set<Index> present_segment_ids;
  for (const Segment* seg : segments) {
    present_segment_ids.insert(seg->getId());
  }

  // SaturateExternalBonds() saturates every atom with hasExternalBond().
  // Clear the external bonds that are internal to this supermolecule
  // beforehand. This only modifies the local QMMolecule copy, not the
  // MD-level data (see QMAtom::clearExternalBond()).
  Index n_original_atoms = qmmol.size();
  for (QMAtom& atom : qmmol) {
    if (!atom.hasExternalBond()) {
      continue;
    }
    if (present_segment_ids.count(atom.getExternalBondPartnerSegmentId()) > 0) {
      atom.clearExternalBond();
    }
  }

  FragmentSaturator::SaturationResult saturation_result =
      FragmentSaturator::SaturateExternalBonds(qmmol);
  QMMolecule relaxed =
      FragmentSaturator::RelaxNewAtoms(saturation_result.mol, n_original_atoms);

  // Fragment atom lists for PODCoupling. Original atoms are classified by
  // id: [0, n_fragment_A_atoms) is fragment A,
  // [n_fragment_A_atoms, n_fragment_AB_atoms) fragment B, linker atoms
  // (mapped last) belong to neither. A saturating H atom inherits the
  // fragment of the atom it saturates
  // (saturation_result.new_atom_parent_ids).
  std::vector<Index> fragment_A_atoms;
  std::vector<Index> fragment_B_atoms;
  for (const QMAtom& atom : relaxed) {
    Index owning_atom_id = atom.getId();
    if (owning_atom_id >= n_original_atoms) {
      owning_atom_id = saturation_result.new_atom_parent_ids[owning_atom_id];
    }
    if (linker_atom_ids.count(owning_atom_id) > 0) {
      continue;
    }
    if (owning_atom_id < n_fragment_A_atoms) {
      fragment_A_atoms.push_back(atom.getId());
    } else if (owning_atom_id < n_fragment_AB_atoms) {
      fragment_B_atoms.push_back(atom.getId());
    }
  }

  Orbitals orbitalsAB;
  orbitalsAB.QMAtoms() = relaxed;

  if (do_dft_input_ || do_dft_run_ || do_dft_parse_) {
    std::string qmpackage_work_dir = work_dir;

    Logger dft_logger(Log::current_level);
    dft_logger.setMultithreading(false);
    dft_logger.setPreface(Log::info, (boost::format("\nDFT INF ...")).str());
    dft_logger.setPreface(Log::error, (boost::format("\nDFT ERR ...")).str());
    dft_logger.setPreface(Log::warning, (boost::format("\nDFT WAR ...")).str());
    dft_logger.setPreface(Log::debug, (boost::format("\nDFT DBG ...")).str());
    std::string package = dftpackage_options_.get("name").as<std::string>();
    std::unique_ptr<QMPackage> qmpackage = QMPackageFactory().Create(package);
    qmpackage->setLog(&dft_logger);
    qmpackage->setRunDir(qmpackage_work_dir);
    qmpackage->Initialize(dftpackage_options_);

    if (do_dft_input_) {
      std::filesystem::create_directories(qmpackage_work_dir);
      // Always use the DFT package's default starting guess. The dimer
      // guess of IQM::EvalJob does not apply: orbitalsAB is the relaxed,
      // H-saturated supermolecule, and the added H atoms have no monomer
      // orbitals. A requested guess is ignored with a warning.
      if (qmpackage->GuessRequested()) {
        XTP_LOG(Log::warning, pLog)
            << "A DFT guess was requested in the dftpackage options, but "
               "IPodCoupling does not support it (the saturating H atoms "
               "have no monomer orbitals). Using the DFT package's default "
               "starting guess instead."
            << std::flush;
      }
      qmpackage->WriteInputFile(orbitalsAB);
    }

    if (do_dft_run_) {
      XTP_LOG(Log::error, pLog) << "Running DFT" << std::flush;
      bool run_dft_status = qmpackage->Run();
      if (!run_dft_status) {
        SetJobToFailed(jres, pLog, qmpackage->getPackageName() + " run failed");
        WriteLoggerToFile(work_dir + "/dft.log", dft_logger);
        return jres;
      }
    }

    if (do_dft_parse_) {
      bool parse_log_status = qmpackage->ParseLogFile(orbitalsAB);
      if (!parse_log_status) {
        SetJobToFailed(jres, pLog, "LOG parsing failed");
        return jres;
      }
      bool parse_orbitals_status = qmpackage->ParseMOsFile(orbitalsAB);
      if (!parse_orbitals_status) {
        SetJobToFailed(jres, pLog, "Orbitals parsing failed");
        return jres;
      }
    }
    qmpackage->CleanUp();
    WriteLoggerToFile(work_dir + "/dft.log", dft_logger);
  } else {
    try {
      orbitalsAB.ReadFromCpt(orbFileAB);
    } catch (std::runtime_error&) {
      SetJobToFailed(jres, pLog,
                     "Do input: failed loading orbitals from " + orbFileAB);
      return jres;
    }
  }

  if (store_dft_) {
    std::filesystem::create_directories(
        std::filesystem::path(orbFileAB).parent_path());
    orbitalsAB.WriteToCpt(orbFileAB);
  }

  // POD coupling on the converged orbitalsAB (parsed above or read from
  // orbFileAB). Guarded by do_podcoupling_: PODCoupling needs converged
  // orbitals with a basis set name, which do not exist after the "input"
  // task alone (an empty name gives a "basis_sets/.xml" error).
  tools::Property job_summary;
  tools::Property& job_output = job_summary.add("output", "");
  if (do_podcoupling_) {
    try {
      // The fragment assignment is otherwise not visible in any output.
      // Printed at Log::debug, i.e. only with -v/--verbose2
      // (application.cc).
      XTP_LOG(Log::debug, pLog)
          << "PODCoupling: fragment_A_atoms (" << fragment_A_atoms.size()
          << " atoms):" << std::flush;
      for (Index a : fragment_A_atoms) {
        XTP_LOG(Log::debug, pLog) << "  " << a << std::flush;
      }
      XTP_LOG(Log::debug, pLog)
          << "PODCoupling: fragment_B_atoms (" << fragment_B_atoms.size()
          << " atoms):" << std::flush;
      for (Index b : fragment_B_atoms) {
        XTP_LOG(Log::debug, pLog) << "  " << b << std::flush;
      }
      if (!linker_atom_ids.empty()) {
        XTP_LOG(Log::debug, pLog)
            << "PODCoupling: linker_atom_ids, excluded from both fragments ("
            << linker_atom_ids.size() << " atoms):" << std::flush;
        for (Index l : linker_atom_ids) {
          XTP_LOG(Log::debug, pLog) << "  " << l << std::flush;
        }
      }

      PODCoupling pod(orbitalsAB, &pLog, fragment_A_atoms, fragment_B_atoms);
      pod.CalculateCouplings(numberofstatesA_, numberofstatesB_);

      // Output format of DFTcoupling::WriteToProperty (dftcoupling.cc):
      // one coupling element with levelA, levelB and j attributes per
      // level pair. The monomer data of DFTcoupling::Addoutput does not
      // exist for POD2, which has no isolated monomer calculations.
      tools::Property& podcoupling_summary = job_output.add(Identify(), "");
      Index homoA = pod.getFragmentAHomoIndex();
      Index lumoA = pod.getFragmentALumoIndex();
      Index homoB = pod.getFragmentBHomoIndex();
      Index lumoB = pod.getFragmentBLumoIndex();
      // homoA/homoB attributes (lumo = homo + 1), as in
      // DFTcoupling::Addoutput; ReadJobFile uses them to select the hole
      // and electron couplings among the coupling elements below.
      podcoupling_summary.setAttribute("homoA", homoA);
      podcoupling_summary.setAttribute("homoB", homoB);
      for (Index levelA = homoA - numberofstatesA_ + 1;
           levelA <= lumoA + numberofstatesA_ - 1; ++levelA) {
        for (Index levelB = homoB - numberofstatesB_ + 1;
             levelB <= lumoB + numberofstatesB_ - 1; ++levelB) {
          double J_hartree = pod.getCouplingElement(levelA, levelB);
          // getCouplingElement() returns Hartree; "j" is written in eV,
          // as in IQM (GetDFTCouplingFromProp, iqm.cc, converts back).
          double J_ev = J_hartree * tools::conv::hrt2ev;
          tools::Property& coupling = podcoupling_summary.add("coupling", "");
          coupling.setAttribute("levelA", levelA);
          coupling.setAttribute("levelB", levelB);
          coupling.setAttribute("j", (boost::format("%1$1.6e") % J_ev).str());
        }
      }
    } catch (std::runtime_error& error) {
      SetJobToFailed(jres, pLog, std::string("PODCoupling: ") + error.what());
      return jres;
    }
  }

  jres.setOutput(job_summary);
  jres.setStatus(Job::COMPLETE);
  return jres;
}

double IPodCoupling::GetPODCouplingFromProp(const tools::Property& podprop,
                                            Index levelA, Index levelB) const {
  for (const tools::Property* state : podprop.Select("coupling")) {
    Index state1 = state->getAttribute<Index>("levelA");
    Index state2 = state->getAttribute<Index>("levelB");
    if (state1 == levelA && state2 == levelB) {
      return state->getAttribute<double>("j") * tools::conv::ev2hrt;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

void IPodCoupling::ReadJobFile(Topology& top) {
  // Structure as in IQM::ReadJobFile (iqm.cc), but with a single
  // podcoupling node and no configurable level map: the HOMO-HOMO
  // coupling is stored as the hole coupling, LUMO-LUMO (lumo = homo + 1)
  // as the electron coupling.
  QMNBList& nblist = top.NBList();
  Index number_of_pairs = nblist.size();
  Index updated_h = 0;
  Index updated_e = 0;
  Index incomplete_jobs = 0;
  Logger log;
  log.setReportLevel(Log::current_level);

  tools::Property xml;
  xml.LoadFromXML(jobfile_);

  for (tools::Property* job : xml.Select("jobs.job")) {
    if (!job->exists("status")) {
      throw std::runtime_error(
          "Jobfile is malformed. <status> tag missing on job.");
    }
    if (job->get("status").as<std::string>() != "COMPLETE" ||
        !job->exists("output")) {
      incomplete_jobs++;
      continue;
    }

    std::vector<Index> id;
    for (tools::Property* segment : job->Select("input.segment")) {
      id.push_back(segment->getAttribute<Index>("id"));
    }
    if (id.size() != 2) {
      throw std::runtime_error(
          "Getting pair ids from jobfile failed, check jobfile.");
    }

    Segment& segA = top.getSegment(id[0]);
    Segment& segB = top.getSegment(id[1]);
    QMPair* qmp = nblist.FindPair(&segA, &segB);
    if (qmp == nullptr) {
      XTP_LOG(Log::error, log)
          << "No pair " << id[0] << ":" << id[1]
          << " found in the neighbor list. Ignoring" << std::flush;
      continue;
    }
    if (qmp->getType() != QMPair::PairType::Hopping) {
      XTP_LOG(Log::error, log) << "WARNING Pair " << qmp->getId()
                               << " is not of any of the "
                                  "Hopping type. Skipping pair"
                               << std::flush;
      continue;
    }

    const tools::Property& pair_property = job->get("output");
    if (!pair_property.exists(Identify())) {
      continue;
    }
    const tools::Property& podprop = pair_property.get(Identify());
    Index homoA = podprop.getAttribute<Index>("homoA");
    Index homoB = podprop.getAttribute<Index>("homoB");
    Index lumoA = homoA + 1;
    Index lumoB = homoB + 1;

    QMStateType hole = QMStateType(QMStateType::Hole);
    double J_hole = GetPODCouplingFromProp(podprop, homoA, homoB);
    if (!std::isnan(J_hole)) {
      qmp->setJeff(J_hole, hole);
      qmp->setJeff2(J_hole * J_hole, hole);
      updated_h++;
    }

    QMStateType electron = QMStateType(QMStateType::Electron);
    double J_electron = GetPODCouplingFromProp(podprop, lumoA, lumoB);
    if (!std::isnan(J_electron)) {
      qmp->setJeff(J_electron, electron);
      qmp->setJeff2(J_electron * J_electron, electron);
      updated_e++;
    }
  }
  XTP_LOG(Log::error, log) << "Pairs [total:updated(e,h)] " << number_of_pairs
                           << ":(" << updated_e << "," << updated_h
                           << ") Incomplete jobs: " << incomplete_jobs << "\n"
                           << std::flush;
  std::cout << log;
}

void IPodCoupling::SetJobToFailed(Job::JobResult& jres, Logger& pLog,
                                  const std::string& errormessage) {
  XTP_LOG(Log::error, pLog) << errormessage << std::flush;
  std::cout << pLog;
  jres.setError(errormessage);
  jres.setStatus(Job::FAILED);
}

void IPodCoupling::WriteLoggerToFile(const std::string& logfile,
                                     Logger& logger) {
  std::ofstream ofs;
  ofs.open(logfile, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + logfile);
  }
  ofs << logger << std::endl;
  ofs.close();
}

}  // namespace xtp
}  // namespace votca
