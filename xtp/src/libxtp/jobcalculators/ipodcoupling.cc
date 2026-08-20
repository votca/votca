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

  // job tasks -- same, established "tasks" pattern as IQM's own
  // (ParseSpecificOptions, iqm.cc), matching its own option names
  // deliberately, for user familiarity, plus one new task specific to
  // this calculator's own, additional purpose (podcoupling).
  std::string tasks_string = options.get(".tasks").as<std::string>();

  // We split either on a space or a comma
  tools::Tokenizer tokenizedTasks(tasks_string, " ,");
  std::vector<std::string> tasks = tokenizedTasks.ToVector();

  do_dft_input_ = std::find(tasks.begin(), tasks.end(), "input") != tasks.end();
  do_dft_run_ = std::find(tasks.begin(), tasks.end(), "dft") != tasks.end();
  do_dft_parse_ = std::find(tasks.begin(), tasks.end(), "parse") != tasks.end();
  do_podcoupling_ =
      std::find(tasks.begin(), tasks.end(), "podcoupling") != tasks.end();

  store_dft_ = options.ifExistsReturnElseReturnDefault<bool>(".store_dft",
                                                             store_dft_);
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
  // Reuses IQM::WriteJobFile's own, exact "one job per neighbor-list
  // pair" logic directly (iqm.cc) -- this calculator's own job
  // generation is identical at this level: which pairs to evaluate is
  // entirely a property of the neighbor list itself, independent of
  // what EvalJob actually does with each pair.
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
      " has no atom with a real, direct external bond toward segment " +
      std::to_string(target_segment_id) +
      " -- FindLinkingSegments() reported a bond path that does not "
      "actually exist at the atom level. This should not be possible in "
      "practice -- check the underlying mapping/checkpoint data.");
}

std::vector<Segment> IPodCoupling::PositionLinkersAlongChain(
    const Topology& top, const Segment& seg1_positioned,
    const std::vector<const Segment*>& linkers,
    const Segment& seg2_positioned) const {
  std::vector<Segment> results;
  if (linkers.empty()) {
    return results;
  }
  // Reserved up front, to its own, exact final size -- genuinely
  // necessary, not just an optimization: previous_positioned, below,
  // holds a direct pointer into results itself, and a
  // std::vector::push_back() that triggers a reallocation would
  // silently invalidate it on the very next iteration.
  results.reserve(linkers.size());

  // Real, direct, one-covalent-bond-at-a-time walk (design worked
  // through directly with the user) -- deliberately NOT any single,
  // whole-chain shift, since each individual linker segment could,
  // in principle, be wrapped to a genuinely different periodic image
  // than its own neighbors in the chain; only walking hop by hop,
  // each one short (a real bond length, far smaller than any
  // realistic PBC box), guarantees the minimum-image convention is
  // unambiguous at every single step.
  const Segment* previous_positioned = &seg1_positioned;
  for (const Segment* current_raw : linkers) {
    const Atom& prev_boundary = FindBoundaryAtomTowardSegment(
        *previous_positioned, current_raw->getId());
    const Atom& curr_boundary = FindBoundaryAtomTowardSegment(
        *current_raw, previous_positioned->getId());

    // Same, real bond-length-based expected-partner-position
    // calculation FragmentSaturator::SaturateExternalBonds already
    // uses (fragmentsaturator.cc) -- where, physically, the current
    // segment's own boundary atom is actually expected to be,
    // according to the previous, already-positioned segment's own
    // side of this same, real bond.
    Eigen::Vector3d expected_partner_pos =
        prev_boundary.getPos() +
        FragmentSaturator::kDefaultCHBondLengthAngstrom *
            tools::conv::ang2bohr * prev_boundary.getExternalBondDirection();

    // The real, direct shift to apply to the WHOLE current segment,
    // so its own boundary atom ends up at (the closest periodic
    // image of) expected_partner_pos -- confirmed directly, by
    // reading BCShortestConnection's own actual implementation
    // (orthorhombicbox.cc: r_j - r_i) before writing this, that
    // PbShortestConnect(r_i, r_j) returns the shift FROM r_i TO r_j,
    // so r_i here must be the atom's own, current position
    // (curr_boundary), and r_j the position it needs to end up at
    // (expected_partner_pos) -- not the other way around, which would
    // silently give the exact opposite, negated shift instead.
    Eigen::Vector3d shift = top.PbShortestConnect(curr_boundary.getPos(),
                                                  expected_partner_pos);

    Segment current_shifted = *current_raw;
    current_shifted.Translate(shift);
    results.push_back(current_shifted);
    previous_positioned = &results.back();
  }

  // Validates the very last hop too, into seg2_positioned -- if
  // FindLinkingSegments() reported a real chain but the very last
  // linker turns out not to actually be bonded to seg2 at all, this
  // throws too, for the exact same reason as every other hop (see
  // FindBoundaryAtomTowardSegment's own header comment). No shift is
  // applied to seg2_positioned itself here -- it is already,
  // separately, correctly positioned (typically via
  // QMPair::Seg2PbCopy(), the existing, established mechanism for
  // the main pair itself).
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

  // Get the pair's own two segment ids from the job -- same,
  // established pattern as IQM::EvalJob (iqm.cc), reused directly,
  // including the optional "qm_geometry" attribute (defaulting to the
  // ground state, "n", when absent).
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
                       std::to_string(ID_B) +
                       " found in the neighbor list.");
    return jres;
  }

  // Same, established path-naming pattern as IQM::EvalJob (iqm.cc),
  // reused directly -- but without orbFileA/orbFileB or eqm_work_dir
  // at all, since those only ever existed to support the dimer-guess
  // mechanism, deliberately skipped entirely here (design confirmed
  // directly with the user: genuinely not useful for this
  // calculator's own use case, where individual fragment/monomer
  // calculations are not wanted or not possible at all). Uses
  // "pairs_ipodcoupling", not IQM's own "pairs_iqm", so the two
  // calculators' own output files never collide if both are run on
  // the same morphology.
  std::string pair_dir =
      (boost::format("%1%%2%%3%%4%%5%") % "pair" % "_" % ID_A % "_" % ID_B)
          .str();
  std::filesystem::path arg_path;
  std::string orbFileAB =
      (arg_path / ipodcoupling_work_dir / "pairs_ipodcoupling" / frame_dir /
       (boost::format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B %
        ".orb")
           .str())
          .generic_string();
  std::string package_append = "workdir_" + Identify();
  std::string work_dir =
      (arg_path / ipodcoupling_work_dir / package_append / frame_dir /
       pair_dir)
          .generic_string();

  // Real, direct, PBC-correct positioning for the two, real pair
  // segments -- reuses QMPair::Seg2PbCopy() directly (qmpair.cc),
  // matching the same, established, PBC-correct path IQM::EvalJob
  // itself already uses when no linkers are involved.
  const Segment* seg1 = pair->Seg1();
  Segment seg2 = pair->Seg2PbCopy();

  // LINKER SEGMENTS -- real, direct discovery + PBC-correct
  // positioning, per the design worked through directly with the
  // user (Topology::FindLinkingSegments/IPodCoupling::
  // PositionLinkersAlongChain, both already built and confirmed
  // compiling earlier this session). Declared here, BEFORE segments
  // itself below, since segments will hold direct pointers into
  // positioned_linkers -- positioned_linkers itself must genuinely
  // outlive every use of segments for the rest of this function.
  std::vector<Segment> positioned_linkers;
  if (include_linkers_) {
    std::vector<const Segment*> linkers =
        top.FindLinkingSegments(*seg1, seg_B);
    if (!linkers.empty()) {
      positioned_linkers =
          PositionLinkersAlongChain(top, *seg1, linkers, seg2);
    }
  }

  std::vector<const Segment*> segments = {seg1, &seg2};
  for (const Segment& linker : positioned_linkers) {
    segments.push_back(&linker);
  }

  QMMolecule qmmol = mapper.map(*seg1, stateA);
  // Real, direct mapped atom count for fragment A -- captured here
  // directly, rather than re-derived later from seg1->size() (the
  // MD-level segment's own atom count), since nothing guarantees
  // these are equal: SegmentMapper::map()'s own size check only
  // confirms the mapping file's own expected atom count matches the
  // mapped result, not that every MD-level atom is mapped at all.
  Index n_fragment_A_atoms = qmmol.size();
  qmmol.AddContainer(mapper.map(seg2, stateB));
  // Real, direct mapped atom count for fragment A + fragment B
  // together, WITHOUT any linker atoms at all -- genuinely different
  // from n_original_atoms below (which does include linker atoms,
  // once mapped) -- needed to correctly bound fragment B's own atom
  // range once linker atoms are also present, since "everything at
  // or past fragment A's own atoms" is no longer synonymous with
  // "fragment B" at all once linkers exist too.
  Index n_fragment_AB_atoms = qmmol.size();

  // Linker segments are mapped at the ground state ("n") always --
  // deliberately no per-linker-segment-type state map at all (unlike
  // IQM's own linker_names, which pairs each linker TYPE with its own
  // QMState) -- this graph-based design needs no such per-type
  // configuration at all, since inclusion itself is already,
  // entirely determined by real bond connectivity alone, not
  // segment type. The ground state is the most physically sensible
  // default for a genuinely neutral, bridging unit not itself
  // directly involved in the actual charge transfer.
  //
  // linker_atom_ids tracks every atom id (within the ORIGINAL,
  // pre-saturation qmmol -- i.e. the same "owning atom id" space
  // fragment_A_atoms/fragment_B_atoms are computed in, below) that
  // belongs to any linker segment -- per direct agreement with the
  // user, linker atoms belong to NEITHER fragment_A_atoms NOR
  // fragment_B_atoms at all.
  std::set<Index> linker_atom_ids;
  for (const Segment& linker : positioned_linkers) {
    Index linker_start_atom_id = qmmol.size();
    qmmol.AddContainer(mapper.map(linker, QMState("n")));
    for (Index i = linker_start_atom_id; i < qmmol.size(); i++) {
      linker_atom_ids.insert(i);
    }
  }

  // Which of this supermolecule's own segment ids are already
  // present -- an external bond whose own partner segment is among
  // these is already satisfied within the supermolecule itself (both
  // sides are already present together), and must NOT be saturated
  // with a new H at all -- this is exactly the check
  // getExternalBondPartnerSegmentId() was built for, earlier this
  // session, together with the user.
  std::set<Index> present_segment_ids;
  for (const Segment* seg : segments) {
    present_segment_ids.insert(seg->getId());
  }

  // FragmentSaturator::SaturateExternalBonds() itself, unconditionally,
  // saturates every hasExternalBond()==true atom in whatever molecule
  // it is given -- it has no notion of "this one is already satisfied,
  // skip it" at all. Reconciled here, on this call's own local copy
  // only (never the underlying, persisted MD-level data -- see
  // QMAtom::clearExternalBond()'s own header comment for exactly why
  // this is safe, worked through directly with the user before
  // implementing this): any atom whose own external-bond partner
  // segment already turns out to be present in this specific
  // supermolecule has its own external bond cleared directly, before
  // SaturateExternalBonds() ever sees it -- at that point, the bond
  // genuinely is no longer external to THIS supermolecule at all.
  Index n_original_atoms = qmmol.size();
  for (QMAtom& atom : qmmol) {
    if (!atom.hasExternalBond()) {
      continue;
    }
    if (present_segment_ids.count(atom.getExternalBondPartnerSegmentId()) >
        0) {
      atom.clearExternalBond();
    }
  }

  FragmentSaturator::SaturationResult saturation_result =
      FragmentSaturator::SaturateExternalBonds(qmmol);
  QMMolecule relaxed = FragmentSaturator::RelaxNewAtoms(
      saturation_result.mol, n_original_atoms);

  // fragment_A_atoms/fragment_B_atoms for PODCoupling's own
  // constructor -- every ORIGINAL atom's own fragment is already
  // known directly from its own id (seg1's own atoms were mapped
  // first, so ids 0..n_fragment_A_atoms-1 are fragment A;
  // n_fragment_A_atoms..n_fragment_AB_atoms-1 are fragment B; any
  // linker atoms, mapped after those two, are excluded from both
  // entirely -- per direct agreement with the user, they belong to
  // NEITHER fragment at all). Every new, saturating H atom must
  // inherit the SAME classification as whichever original atom it is
  // saturating -- not recoverable from its own final position/index
  // alone (a real gap surfaced and worked through directly with the
  // user before extending SaturateExternalBonds's own return type to
  // track this explicitly, saturation_result.new_atom_parent_ids).
  std::vector<Index> fragment_A_atoms;
  std::vector<Index> fragment_B_atoms;
  for (const QMAtom& atom : relaxed) {
    Index owning_atom_id = atom.getId();
    if (owning_atom_id >= n_original_atoms) {
      owning_atom_id =
          saturation_result.new_atom_parent_ids[owning_atom_id];
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
    dft_logger.setPreface(Log::warning,
                          (boost::format("\nDFT WAR ...")).str());
    dft_logger.setPreface(Log::debug, (boost::format("\nDFT DBG ...")).str());
    std::string package = dftpackage_options_.get("name").as<std::string>();
    std::unique_ptr<QMPackage> qmpackage = QMPackageFactory().Create(package);
    qmpackage->setLog(&dft_logger);
    qmpackage->setRunDir(qmpackage_work_dir);
    qmpackage->Initialize(dftpackage_options_);

    if (do_dft_input_) {
      std::filesystem::create_directories(qmpackage_work_dir);
      // Deliberately, always the plain "no guess, start from the DFT
      // package's own default starting guess" path -- the dimer-guess
      // mechanism IQM::EvalJob itself optionally uses (combining two,
      // separately pre-computed monomer orbital files) is skipped
      // entirely here, on purpose (design confirmed directly with the
      // user): our own orbitalsAB is not a simple two-monomer
      // combination at all -- it is the RELAXED, H-SATURATED
      // supermolecule, and the new H atom(s) have no corresponding
      // monomer orbitals to guess from in the first place. If the DFT
      // package's own options actually requested a guess anyway, warn
      // directly rather than silently ignoring the request.
      if (qmpackage->GuessRequested()) {
        XTP_LOG(Log::warning, pLog)
            << "A DFT guess was requested in the dftpackage options, but "
               "IPodCoupling does not support this at all (its own "
               "supermolecule is not a simple monomer combination, given "
               "the new, saturating H atom(s) have no corresponding "
               "monomer orbitals to guess from) -- proceeding with the "
               "DFT package's own default starting guess instead."
            << std::flush;
      }
      qmpackage->WriteInputFile(orbitalsAB);
    }

    if (do_dft_run_) {
      XTP_LOG(Log::error, pLog) << "Running DFT" << std::flush;
      bool run_dft_status = qmpackage->Run();
      if (!run_dft_status) {
        SetJobToFailed(jres, pLog,
                       qmpackage->getPackageName() + " run failed");
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

  // Real PODCoupling calculation, on the real, converged orbitalsAB
  // (either freshly computed above, via do_dft_parse_, or read back
  // from a previously-stored orbFileAB) -- final piece of the
  // six-step design worked through directly with the user at the
  // very start of this whole calculator.
  //
  // Real, direct bug fix: this whole block used to run
  // unconditionally, with no do_podcoupling_ guard at all -- caught
  // directly by the user's own real, direct run (tasks="input" only,
  // to inspect the DFT input file directly): PODCoupling's own
  // constructor genuinely needs a real, converged orbitalsAB with an
  // actual basis set name set on it, which simply does not exist yet
  // at all when only the DFT input file has been written (do_dft_run_/
  // do_dft_parse_ both false) -- attempting it anyway threw a real,
  // confusing "basis_sets/.xml" error (an empty basis set name) even
  // though podcoupling was never requested as a task at all.
  tools::Property job_summary;
  tools::Property& job_output = job_summary.add("output", "");
  if (do_podcoupling_) {
    try {
      // Real, direct debug-level output, worked through directly with
      // the user -- fragment_A_atoms/fragment_B_atoms (PODCoupling's
      // own terminology for the same thing VOTCA itself calls a
      // "segment") are otherwise entirely internal: unlike the DFT
      // input file's own, directly-visible written coordinates, there
      // is no other, direct way for a user to confirm which real atom
      // indices actually ended up in which fragment at all, short of
      // running a real, full, genuinely expensive DFT calculation
      // first (do_podcoupling_ itself already requires a real,
      // converged orbitalsAB with a real basis set name set on it --
      // confirmed directly, earlier this same session, per the
      // do_podcoupling_ guard fix). Printed at Log::debug specifically
      // (not error/warning/info) -- confirmed directly, from
      // application.cc, that this is exactly the level -v/--verbose2
      // itself raises current_level to; the two milder --verbose/
      // --verbose1 flags only ever reach warning/info, not this.
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

      // Same, established per-pair output format as DFTcoupling's own
      // WriteToProperty (dftcoupling.cc), confirmed directly by reading
      // it before writing this, rather than invented separately --
      // <coupling levelA="..." levelB="..." j="..."/> for every
      // (levelA, levelB) pair within the requested range, matching
      // DFTcoupling's own, backward-compatible core format exactly.
      // Unlike DFTcoupling's own, more elaborate Addoutput (monomer
      // energies, raw TB matrices, diagnostics), none of that is added
      // here at all -- it genuinely does not apply to POD2 at all,
      // which has no separate, isolated monomer calculation to compare
      // against in the first place (the whole point of POD2, per this
      // class's own header comment, podcoupling.h).
      tools::Property& podcoupling_summary = job_output.add(Identify(), "");
      Index homoA = pod.getFragmentAHomoIndex();
      Index lumoA = pod.getFragmentALumoIndex();
      Index homoB = pod.getFragmentBHomoIndex();
      Index lumoB = pod.getFragmentBLumoIndex();
      // homoA/homoB (lumoA/lumoB = homoA/homoB + 1, always -- see
      // PODCoupling::getFragmentALumoIndex's own header comment,
      // podcoupling.h) written directly as attributes on this same
      // node, matching DFTcoupling::Addoutput's own, established
      // pattern (dftcoupling.cc: dftcoupling.setAttribute("homoA",
      // orbitalsA.getHomo())) exactly -- genuinely needed for
      // ReadJobFile to be able to translate "the hole/electron
      // coupling" into the actual, specific levelA/levelB pair among
      // the (potentially many, if numberofstatesA_/B_ > 1) <coupling>
      // elements written below.
      podcoupling_summary.setAttribute("homoA", homoA);
      podcoupling_summary.setAttribute("homoB", homoB);
      for (Index levelA = homoA - numberofstatesA_ + 1;
          levelA <= lumoA + numberofstatesA_ - 1; ++levelA) {
        for (Index levelB = homoB - numberofstatesB_ + 1;
            levelB <= lumoB + numberofstatesB_ - 1; ++levelB) {
          double J_hartree = pod.getCouplingElement(levelA, levelB);
          // Written in eV, not Hartree -- matching IQM::WriteToProperty's
          // own, established convention exactly (confirmed directly by
          // reading GetDFTCouplingFromProp, iqm.cc, which converts the
          // read-back "j" value FROM eV back TO Hartree via
          // tools::conv::ev2hrt -- meaning IQM's own "j" is written in
          // eV, even though PODCoupling::getCouplingElement's own
          // documented return unit, podcoupling.h, is Hartree). More
          // human-readable too -- typical couplings are meV-scale, not
          // the much smaller numbers raw Hartree would give.
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
                                            Index levelA,
                                            Index levelB) const {
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
  // Same, established structure as IQM::ReadJobFile (iqm.cc),
  // confirmed directly by reading it in full before writing this,
  // rather than guessed -- but reads back a single, simpler
  // "podcoupling" node instead of IQM's own separate dftcoupling/
  // bsecoupling nodes, and, unlike IQM, has no hole_levels_/
  // electron_levels_-style, user-configurable state map at all: the
  // HOMO-HOMO coupling is always written back as the hole coupling,
  // and LUMO-LUMO as the electron coupling -- the most direct,
  // standard headline result, matching PODCoupling's own
  // getFragmentAHomoIndex()/getFragmentALumoIndex() (podcoupling.h)
  // convention exactly (lumo = homo + 1, always).
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
