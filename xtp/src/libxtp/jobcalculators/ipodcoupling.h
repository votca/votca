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

#pragma once
#ifndef VOTCA_XTP_IPODCOUPLING_H
#define VOTCA_XTP_IPODCOUPLING_H

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/parallelxjobcalc.h"

namespace votca {
namespace xtp {

/**
 * \brief H-saturated, DFT + POD2-based electronic coupling calculator
 *
 * A new, separate jobcalculator (deliberately not a modification of
 * IQM itself -- IQM already has an established, widely-used role, and
 * directly entangling this genuinely different workflow into it would
 * risk breaking existing behavior for everyone who already relies on
 * it, per the design worked through directly with the user before
 * starting this implementation).
 *
 * For each neighbor-list pair: maps both segments (plus, in future,
 * any linking segments -- see this class's own EvalJob header comment
 * for the current, still-open status of that piece), combines them
 * into one supermolecule, saturates any remaining covalently-cut
 * fragment boundaries with new H atoms (FragmentSaturator::
 * SaturateExternalBonds, geometric placement, then RelaxNewAtoms, a
 * constrained OpenBabel MMFF94/UFF relaxation -- both already built
 * and tested earlier this session), runs the actual DFT calculation
 * on the resulting, H-saturated geometry, then runs PODCoupling on
 * the resulting orbitals.
 *
 * Callname: ipodcoupling
 */
class IPodCoupling final : public ParallelXJobCalc<std::vector<Job> > {
 public:
  std::string Identify() const { return "ipodcoupling"; }
  Job::JobResult EvalJob(const Topology& top, Job& job, QMThread& opThread);
  void WriteJobFile(const Topology& top);
  void ReadJobFile(Topology& top);

 protected:
  void ParseSpecificOptions(const tools::Property& user_options);

 private:
  void SetJobToFailed(Job::JobResult& jres, Logger& pLog,
                      const std::string& errormessage);
  void WriteLoggerToFile(const std::string& logfile, Logger& logger);

  // Finds, within podprop (the same "podcoupling" node EvalJob writes
  // via podcoupling_summary), the <coupling> element with the given,
  // exact levelA/levelB, and returns its own "j" value converted back
  // to Hartree (podprop's own "j" is written in eV -- see EvalJob's
  // own comment for exactly why) -- NaN if no such element exists at
  // all (e.g. this specific pair was never actually within the
  // requested numberofstatesA_/B_ range at the time the job was
  // originally run). Same, established pattern as IQM::
  // GetDFTCouplingFromProp (iqm.cc), confirmed directly by reading it
  // before writing this.
  double GetPODCouplingFromProp(const tools::Property& podprop, Index levelA,
                                Index levelB) const;

  // Finds, within seg, the specific atom whose own external bond
  // points toward target_segment_id -- genuinely different from just
  // "any atom with hasExternalBond()", since a segment can, in
  // principle, have external bonds toward more than one other
  // segment at all (only the one specifically bonded toward
  // target_segment_id is wanted here). Throws directly (per direct
  // agreement with the user) if no such atom exists at all -- this
  // would mean seg and the segment target_segment_id refers to are
  // NOT actually, directly bonded, despite FindLinkingSegments()
  // itself already having reported a real bond path between them; a
  // silent skip/fallback here would produce a geometrically wrong
  // supermolecule with no warning at all, worse than failing loudly.
  const Atom& FindBoundaryAtomTowardSegment(const Segment& seg,
                                            Index target_segment_id) const;

  // Returns PBC-correctly-positioned COPIES of every segment in
  // linkers (never the originals stored in Topology, matching
  // QMPair::Seg2PbCopy()'s own, established "return a fresh copy"
  // convention exactly) -- walks the real, actual covalent chain
  // seg1_positioned -> linkers[0] -> linkers[1] -> ... ->
  // seg2_positioned one bond at a time (see this method's own .cc
  // implementation for exactly why one bond at a time, not any
  // single, whole-chain shift, is needed), using
  // Topology::PbShortestConnect() at each hop, anchored on each
  // hop's own, specific FindBoundaryAtomTowardSegment result -- not
  // any segment's own, whole-molecule center of mass, unlike
  // QMPair::Seg2PbCopy() itself, since a large segment's own center
  // can be genuinely far from its own, specific bonded atom.
  //
  // seg1_positioned/seg2_positioned are assumed to already be
  // correctly positioned themselves (this walk's own two, fixed
  // anchors) -- typically seg1 itself (unmoved, the pair's own,
  // established anchor) and pair->Seg2PbCopy() (already
  // PBC-corrected, via the existing, established pair-level
  // mechanism) respectively. The very last hop
  // (linkers.back() -> seg2_positioned) is validated the exact same
  // way as every other hop -- if FindLinkingSegments() reported a
  // real chain but the very last linker turns out not to be actually
  // bonded to seg2 at all, this throws too, for the exact same
  // reason.
  std::vector<Segment> PositionLinkersAlongChain(
      const Topology& top, const Segment& seg1_positioned,
      const std::vector<const Segment*>& linkers,
      const Segment& seg2_positioned) const;

  tools::Property dftpackage_options_;
  tools::Property podcoupling_options_;

  // PODCoupling::CalculateCouplings's own numberofstatesA/B --
  // matches podcoupling.xml's own, established "levA"/"levB" option
  // names exactly (already used by the standalone podcoupling tool,
  // podcouplingtool.cc), for user familiarity -- unlike that tool,
  // fragment_A/fragment_B themselves are never parsed as options at
  // all here; they are computed directly from the segment mapping
  // itself (IPodCoupling::EvalJob), not user-supplied.
  Index numberofstatesA_ = 1;
  Index numberofstatesB_ = 1;

  // Whether to look for, and include, linking segments at all (via
  // Topology::FindLinkingSegments) -- deliberately a simple,
  // real toggle, per direct agreement with the user (default off,
  // matching IQM's own default-off linker behavior). Unlike IQM's
  // own linker_names (a segment-TYPE whitelist), this graph-based
  // design needs no such whitelist at all -- inclusion is already,
  // entirely determined by real bond connectivity alone (is there an
  // actual, real covalent path between seg_A and seg_B at all).
  // Deliberately no cap on the number of linker segments included at
  // all either, per direct agreement with the user.
  bool include_linkers_ = false;

  // What to do -- same, established "tasks" pattern as IQM's own
  // do_dft_input_/do_dft_run_/do_dft_parse_ (matching option names,
  // deliberately, for user familiarity), plus a new task specific to
  // this calculator's own, additional purpose.
  bool do_dft_input_ = false;
  bool do_dft_run_ = false;
  bool do_dft_parse_ = false;
  bool do_podcoupling_ = false;

  bool store_dft_ = false;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_IPODCOUPLING_H
