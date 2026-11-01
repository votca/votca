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
 * Kept separate from IQM so that IQM's existing behavior is unchanged.
 *
 * For each neighbor-list pair: maps both segments (and, optionally, the
 * linking segments between them) into one supermolecule, saturates the
 * remaining covalently cut boundaries with H atoms
 * (FragmentSaturator::SaturateExternalBonds, then RelaxNewAtoms, a
 * constrained OpenBabel MMFF94/UFF relaxation), runs DFT on the saturated
 * geometry and evaluates PODCoupling on the resulting orbitals.
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

  // Returns the "j" value (stored in eV) of the coupling element of podprop
  // with the given levelA/levelB, converted to Hartree; NaN if the levels
  // were outside the range computed in EvalJob. As in
  // IQM::GetDFTCouplingFromProp (iqm.cc).
  double GetPODCouplingFromProp(const tools::Property& podprop, Index levelA,
                                Index levelB) const;

  // Returns the atom of seg whose external bond points to segment
  // target_segment_id (a segment may have external bonds to several
  // segments). Throws if there is none, i.e. if the bond path reported by
  // FindLinkingSegments() does not exist at the atom level; a silent
  // fallback would produce a wrong supermolecule geometry.
  const Atom& FindBoundaryAtomTowardSegment(const Segment& seg,
                                            Index target_segment_id) const;

  // Returns PBC-corrected copies of the linker segments (the Topology is
  // not modified, as for QMPair::Seg2PbCopy()). Walks the covalent chain
  // seg1_positioned -> linkers[0] -> ... -> seg2_positioned one bond at a
  // time with Topology::PbShortestConnect(), anchored on the bonded
  // boundary atoms (FindBoundaryAtomTowardSegment) rather than on centers
  // of mass, which can be far from the bond in large segments.
  //
  // seg1_positioned and seg2_positioned must already be positioned
  // (typically seg1 and pair->Seg2PbCopy()). The last hop to
  // seg2_positioned is validated like all others and throws if broken.
  std::vector<Segment> PositionLinkersAlongChain(
      const Topology& top, const Segment& seg1_positioned,
      const std::vector<const Segment*>& linkers,
      const Segment& seg2_positioned) const;

  tools::Property dftpackage_options_;
  tools::Property podcoupling_options_;

  // numberofstatesA/B of PODCoupling::CalculateCouplings, read from the
  // "levA"/"levB" options as in the podcoupling tool. Unlike that tool,
  // the fragments are derived from the segment mapping, not user input.
  Index numberofstatesA_ = 1;
  Index numberofstatesB_ = 1;

  // Include linking segments found by Topology::FindLinkingSegments
  // (default off, as in IQM). Selection is by covalent connectivity
  // between seg_A and seg_B, not by a segment-type list as in IQM's
  // linker_names; the number of linkers is not capped.
  bool include_linkers_ = false;

  // Tasks, as in IQM, plus "podcoupling".
  bool do_dft_input_ = false;
  bool do_dft_run_ = false;
  bool do_dft_parse_ = false;
  bool do_podcoupling_ = false;

  bool store_dft_ = false;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_IPODCOUPLING_H
