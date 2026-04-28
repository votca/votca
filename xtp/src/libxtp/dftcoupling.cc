/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/dftcoupling.h"

namespace votca {
namespace xtp {

using boost::format;
using std::flush;

void DFTcoupling::Initialize(tools::Property& options) {
  degeneracy_ = options.get("degeneracy").as<double>();
  degeneracy_ *= tools::conv::ev2hrt;
  numberofstatesA_ = options.get("levA").as<Index>();
  numberofstatesB_ = options.get("levB").as<Index>();
  // Optional TB output: monomer MO energies (KS and QP), raw H/S matrices,
  // and diagnostics. Default false to keep XML output compact for KMC/rate
  // workflows that only need the scalar effective couplings.
  output_tb_ =
      options.ifExistsReturnElseReturnDefault<bool>("output_tb", false);
}

// =============================================================================
// WriteMatrixToProperty
// Write an Eigen matrix as space-separated rows into an XML property node.
// Each row is stored as attribute "row_N". Unit conversion applied if given.
// =============================================================================
void DFTcoupling::WriteMatrixToProperty(tools::Property& prop,
                                        const std::string& name,
                                        const Eigen::MatrixXd& mat,
                                        double conversion) {
  tools::Property& mat_prop = prop.add(name, "");
  mat_prop.setAttribute("rows", mat.rows());
  mat_prop.setAttribute("cols", mat.cols());
  for (Index i = 0; i < mat.rows(); i++) {
    std::string row_str;
    for (Index j = 0; j < mat.cols(); j++) {
      if (j > 0) row_str += " ";
      row_str += (format("%1$1.6e") % (mat(i, j) * conversion)).str();
    }
    mat_prop.setAttribute("row_" + std::to_string(i), row_str);
  }
}

// =============================================================================
// WriteToProperty
// Write the scalar effective coupling for one (levelA, levelB) pair.
// Backward compatible: attribute names unchanged.
// =============================================================================
void DFTcoupling::WriteToProperty(tools::Property& type_summary,
                                  const Orbitals& orbitalsA,
                                  const Orbitals& orbitalsB, Index a,
                                  Index b) const {
  double J = getCouplingElement(a, b, orbitalsA, orbitalsB);
  tools::Property& coupling = type_summary.add("coupling", "");
  coupling.setAttribute("levelA", a);
  coupling.setAttribute("levelB", b);
  coupling.setAttribute("j", (format("%1$1.6e") % J).str());
}

// =============================================================================
// Addoutput
// Write all results to the XML property tree.
// Extends the existing pairwise coupling output with:
//   - monomer_energies: isolated MO energies [eV] for TB site energies
//   - tb_matrices: raw H and S block matrices [H in eV, S dimensionless]
//   - diagnostics: minimum S eigenvalue as basis quality indicator
// The existing <coupling> elements are unchanged for backward compatibility.
// =============================================================================
void DFTcoupling::Addoutput(tools::Property& type_summary,
                            const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB) const {
  tools::Property& dftcoupling = type_summary.add(Identify(), "");
  dftcoupling.setAttribute("homoA", orbitalsA.getHomo());
  dftcoupling.setAttribute("homoB", orbitalsB.getHomo());

  // helper lambda: write one carrier type (hole or electron) into a
  // property node, including existing pairwise couplings, monomer energies
  // (KS always, QPpert when available), raw TB matrices, and diagnostics.
  auto WriteCarrier = [&](tools::Property& carrier_summary, Index start_a,
                          Index end_a, Index start_b, Index end_b,
                          const Eigen::MatrixXd& JAB_dimer,
                          const Eigen::MatrixXd& S_AxB,
                          const Eigen::VectorXd& moEnA_KS,
                          const Eigen::VectorXd& moEnB_KS,
                          const Eigen::VectorXd& moEnA_QP,
                          const Eigen::VectorXd& moEnB_QP, double min_S_eval) {
    // --- existing pairwise effective couplings (backward compatible) ---
    for (Index a = start_a; a <= end_a; ++a) {
      for (Index b = start_b; b <= end_b; ++b) {
        WriteToProperty(carrier_summary, orbitalsA, orbitalsB, a, b);
      }
    }

    // --- TB output: monomer energies, raw matrices, diagnostics ---
    // Only written when output_tb_ = true. For KMC/rate workflows that only
    // need scalar couplings, leave output_tb_ = false to keep XML compact.
    if (output_tb_) {
      // Monomer MO energies
      // eV_KS: Kohn-Sham DFT eigenvalue (always present)
      // eV_QP: QPpert quasiparticle energy (only when GW was run)
      {
        bool hasQP_A = (moEnA_QP.size() == moEnA_KS.size());
        bool hasQP_B = (moEnB_QP.size() == moEnB_KS.size());

        tools::Property& mono_prop =
            carrier_summary.add("monomer_energies", "");
        mono_prop.setAttribute("has_qp",
                               (hasQP_A && hasQP_B) ? "true" : "false");

        tools::Property& propA = mono_prop.add("fragmentA", "");
        for (Index i = 0; i < moEnA_KS.size(); i++) {
          tools::Property& e = propA.add("mo", "");
          e.setAttribute("level", start_a + i);
          e.setAttribute("eV_KS", (format("%1$1.6e") % moEnA_KS(i)).str());
          if (hasQP_A) {
            e.setAttribute("eV_QP", (format("%1$1.6e") % moEnA_QP(i)).str());
          }
        }
        tools::Property& propB = mono_prop.add("fragmentB", "");
        for (Index i = 0; i < moEnB_KS.size(); i++) {
          tools::Property& e = propB.add("mo", "");
          e.setAttribute("level", start_b + i);
          e.setAttribute("eV_KS", (format("%1$1.6e") % moEnB_KS(i)).str());
          if (hasQP_B) {
            e.setAttribute("eV_QP", (format("%1$1.6e") % moEnB_QP(i)).str());
          }
        }
      }

      // Raw TB matrices (H in eV, S dimensionless)
      {
        Index levA_range = moEnA_KS.size();
        Index levB_range = moEnB_KS.size();
        tools::Property& tb_prop = carrier_summary.add("tb_matrices", "");
        tb_prop.setAttribute("n_A", levA_range);
        tb_prop.setAttribute("n_B", levB_range);

        WriteMatrixToProperty(tb_prop, "H_AA",
                              JAB_dimer.topLeftCorner(levA_range, levA_range),
                              votca::tools::conv::hrt2ev);
        WriteMatrixToProperty(tb_prop, "S_AA",
                              S_AxB.topLeftCorner(levA_range, levA_range));
        WriteMatrixToProperty(
            tb_prop, "H_BB",
            JAB_dimer.bottomRightCorner(levB_range, levB_range),
            votca::tools::conv::hrt2ev);
        WriteMatrixToProperty(tb_prop, "S_BB",
                              S_AxB.bottomRightCorner(levB_range, levB_range));
        WriteMatrixToProperty(tb_prop, "H_AB",
                              JAB_dimer.topRightCorner(levA_range, levB_range),
                              votca::tools::conv::hrt2ev);
        WriteMatrixToProperty(tb_prop, "S_AB",
                              S_AxB.topRightCorner(levA_range, levB_range));
      }

      // Basis diagnostics
      {
        tools::Property& diag_prop = carrier_summary.add("diagnostics", "");
        diag_prop.setAttribute("min_S_eigenvalue",
                               (format("%1$1.6e") % min_S_eval).str());
        diag_prop.setAttribute("basis_ok",
                               (min_S_eval > 0.01) ? "true" : "false");
      }
    }  // output_tb_
  };

  // --- hole block ---
  tools::Property& hole_summary = dftcoupling.add("hole", "");
  WriteCarrier(hole_summary, Range_orbA.first, orbitalsA.getHomo(),
               Range_orbB.first, orbitalsB.getHomo(), JAB_dimer_hole_,
               S_AxB_hole_, moEnergiesA_hole_KS_, moEnergiesB_hole_KS_,
               moEnergiesA_hole_QP_, moEnergiesB_hole_QP_,
               min_S_eigenvalue_hole_);

  // --- electron block ---
  tools::Property& electron_summary = dftcoupling.add("electron", "");
  WriteCarrier(electron_summary, orbitalsA.getLumo(),
               Range_orbA.first + Range_orbA.second - 1, orbitalsB.getLumo(),
               Range_orbB.first + Range_orbB.second - 1, JAB_dimer_elec_,
               S_AxB_elec_, moEnergiesA_elec_KS_, moEnergiesB_elec_KS_,
               moEnergiesA_elec_QP_, moEnergiesB_elec_QP_,
               min_S_eigenvalue_elec_);
}

std::pair<Index, Index> DFTcoupling::DetermineRangeOfStates(
    const Orbitals& orbital, Index numberofstates) const {
  const Eigen::VectorXd& MOEnergies = orbital.MOs().eigenvalues();
  if (std::abs(MOEnergies(orbital.getHomo()) - MOEnergies(orbital.getLumo())) <
      degeneracy_) {
    throw std::runtime_error(
        "Homo Lumo Gap is smaller than degeneracy. "
        "Either your degeneracy is too large or your Homo and Lumo are "
        "degenerate");
  }

  Index minimal = orbital.getHomo() - numberofstates + 1;
  Index maximal = orbital.getLumo() + numberofstates - 1;

  std::vector<Index> deg_min = orbital.CheckDegeneracy(minimal, degeneracy_);
  minimal = *std::min_element(deg_min.begin(), deg_min.end());

  std::vector<Index> deg_max = orbital.CheckDegeneracy(maximal, degeneracy_);
  maximal = *std::max_element(deg_max.begin(), deg_max.end());

  std::pair<Index, Index> result;
  result.first = minimal;                 // start index
  result.second = maximal - minimal + 1;  // size

  return result;
}

double DFTcoupling::getCouplingElement(Index levelA, Index levelB,
                                       const Orbitals& orbitalsA,
                                       const Orbitals& orbitalsB) const {
  Index levelsA = Range_orbA.second;
  if (degeneracy_ != 0) {
    std::vector<Index> list_levelsA =
        orbitalsA.CheckDegeneracy(levelA, degeneracy_);
    std::vector<Index> list_levelsB =
        orbitalsB.CheckDegeneracy(levelB, degeneracy_);

    double JAB_sq = 0;
    for (Index iA : list_levelsA) {
      Index indexA = iA - Range_orbA.first;
      for (Index iB : list_levelsB) {
        Index indexB = iB - Range_orbB.first + levelsA;
        double JAB_one_level = JAB(indexA, indexB);
        JAB_sq += JAB_one_level * JAB_one_level;
      }
    }
    return std::sqrt(JAB_sq /
                     double(list_levelsA.size() * list_levelsB.size())) *
           tools::conv::hrt2ev;
  } else {
    Index indexA = levelA - Range_orbA.first;
    Index indexB = levelB - Range_orbB.first + levelsA;
    return JAB(indexA, indexB) * tools::conv::hrt2ev;
  }
}

// =============================================================================
// CalculateCouplings
// Computes the Lowdin-orthogonalized effective coupling matrix JAB, and also
// stores the raw (pre-Lowdin) projected Fock matrix and overlap for each
// carrier type, along with isolated monomer MO energies, for TB use.
// =============================================================================
void DFTcoupling::CalculateCouplings(const Orbitals& orbitalsA,
                                     const Orbitals& orbitalsB,
                                     const Orbitals& orbitalsAB) {

  XTP_LOG(Log::error, *pLog_) << "Calculating electronic couplings" << flush;

  CheckAtomCoordinates(orbitalsA, orbitalsB, orbitalsAB);

  Index basisA = orbitalsA.getBasisSetSize();
  Index basisB = orbitalsB.getBasisSetSize();

  if ((basisA == 0) || (basisB == 0)) {
    throw std::runtime_error("Basis set size is not stored in monomers");
  }

  Range_orbA = DetermineRangeOfStates(orbitalsA, numberofstatesA_);
  Range_orbB = DetermineRangeOfStates(orbitalsB, numberofstatesB_);

  Index levelsA = Range_orbA.second;
  Index levelsB = Range_orbB.second;

  XTP_LOG(Log::error, *pLog_)
      << "Levels:Basis A[" << levelsA << ":" << basisA << "]"
      << " B[" << levelsB << ":" << basisB << "]" << flush;

  if ((levelsA == 0) || (levelsB == 0)) {
    throw std::runtime_error(
        "No information about number of occupied/unoccupied levels is stored");
  }

  // --- project monomer MOs onto dimer orbital basis ---
  auto MOsA = orbitalsA.MOs().eigenvectors().middleCols(Range_orbA.first,
                                                        Range_orbA.second);
  auto MOsB = orbitalsB.MOs().eigenvectors().middleCols(Range_orbB.first,
                                                        Range_orbB.second);

  XTP_LOG(Log::info, *pLog_) << "Calculating overlap matrix for basisset: "
                             << orbitalsAB.getDFTbasisName() << flush;

  Eigen::MatrixXd overlap =
      CalculateOverlapMatrix(orbitalsAB) * orbitalsAB.MOs().eigenvectors();

  XTP_LOG(Log::info, *pLog_)
      << "Projecting monomers onto dimer orbitals" << flush;
  Eigen::MatrixXd A_AB = MOsA.transpose() * overlap.topRows(basisA);
  Eigen::MatrixXd B_AB = MOsB.transpose() * overlap.bottomRows(basisB);

  Eigen::VectorXd mag_A = A_AB.rowwise().squaredNorm();
  if (mag_A.any() < 0.95) {
    XTP_LOG(Log::error, *pLog_)
        << "\nWarning: "
        << "Projection of orbitals of monomer A on dimer is insufficient,mag="
        << mag_A.minCoeff() << flush;
  }
  Eigen::VectorXd mag_B = B_AB.rowwise().squaredNorm();
  if (mag_B.any() < 0.95) {
    XTP_LOG(Log::error, *pLog_)
        << "\nWarning: "
        << "Projection of orbitals of monomer B on dimer is insufficient,mag="
        << mag_B.minCoeff() << flush;
  }

  // --- stack projected MOs into combined basis ---
  Eigen::MatrixXd psi_AxB_dimer_basis(levelsA + levelsB, A_AB.cols());
  psi_AxB_dimer_basis.topRows(levelsA) = A_AB;
  psi_AxB_dimer_basis.bottomRows(levelsB) = B_AB;

  // --- full projected Fock matrix and overlap (all MOs combined) ---
  XTP_LOG(Log::info, *pLog_)
      << "Projecting the Fock matrix onto the dimer basis" << flush;
  Eigen::MatrixXd JAB_dimer_full = psi_AxB_dimer_basis *
                                   orbitalsAB.MOs().eigenvalues().asDiagonal() *
                                   psi_AxB_dimer_basis.transpose();

  XTP_LOG(Log::info, *pLog_) << "Constructing Overlap matrix" << flush;
  Eigen::MatrixXd S_AxB_full =
      psi_AxB_dimer_basis * psi_AxB_dimer_basis.transpose();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_AxB_full);
  Eigen::MatrixXd Sm1 = es.operatorInverseSqrt();
  XTP_LOG(Log::info, *pLog_) << "Smallest eigenvalue of overlap matrix is "
                             << es.eigenvalues()(0) << flush;

  // Lowdin-orthogonalized effective coupling (existing, unchanged)
  JAB = Sm1 * JAB_dimer_full * Sm1;

  // -----------------------------------------------------------------------
  // New: store raw matrices and monomer energies per carrier type.
  //
  // The combined basis is ordered [A_occ, ..., A_unocc, B_occ, ..., B_unocc]
  // since Range_orbA and Range_orbB each span HOMO-n to LUMO+n.
  //
  // For TB purposes we want separate hole and electron sub-blocks.
  // Within the full (levelsA + levelsB) basis:
  //   hole A:     rows/cols  [0, holesA)
  //   electron A: rows/cols  [holesA, levelsA)
  //   hole B:     rows/cols  [levelsA, levelsA+holesB)
  //   electron B: rows/cols  [levelsA+holesB, levelsA+levelsB)
  // where holesA = homoA - Range_orbA.first + 1, etc.
  // -----------------------------------------------------------------------

  Index homoA_idx = orbitalsA.getHomo();
  Index homoB_idx = orbitalsB.getHomo();

  // Number of hole and electron states per fragment within the range
  Index holesA = homoA_idx - Range_orbA.first + 1;
  Index elecsA = levelsA - holesA;
  Index holesB = homoB_idx - Range_orbB.first + 1;
  Index elecsB = levelsB - holesB;

  // Permutation to reorder [A_occ, A_unocc, B_occ, B_unocc] into
  // [A_occ, B_occ] for holes and [A_unocc, B_unocc] for electrons.
  // Build index vectors for each block.
  std::vector<Index> hole_idx, elec_idx;
  for (Index i = 0; i < holesA; i++) hole_idx.push_back(i);
  for (Index i = 0; i < holesB; i++) hole_idx.push_back(levelsA + i);
  for (Index i = holesA; i < levelsA; i++) elec_idx.push_back(i);
  for (Index i = holesB; i < levelsB; i++) elec_idx.push_back(levelsA + i);

  // Extract hole sub-block
  Index nh = hole_idx.size();
  Eigen::MatrixXd JAB_hole(nh, nh), S_hole(nh, nh);
  for (Index i = 0; i < nh; i++) {
    for (Index j = 0; j < nh; j++) {
      JAB_hole(i, j) = JAB_dimer_full(hole_idx[i], hole_idx[j]);
      S_hole(i, j) = S_AxB_full(hole_idx[i], hole_idx[j]);
    }
  }

  // Extract electron sub-block
  Index ne = elec_idx.size();
  Eigen::MatrixXd JAB_elec(ne, ne), S_elec(ne, ne);
  for (Index i = 0; i < ne; i++) {
    for (Index j = 0; j < ne; j++) {
      JAB_elec(i, j) = JAB_dimer_full(elec_idx[i], elec_idx[j]);
      S_elec(i, j) = S_AxB_full(elec_idx[i], elec_idx[j]);
    }
  }

  // Store raw matrices
  JAB_dimer_hole_ = JAB_hole;
  S_AxB_hole_ = S_hole;
  JAB_dimer_elec_ = JAB_elec;
  S_AxB_elec_ = S_elec;

  // Minimum S eigenvalues per carrier type
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_h(S_hole);
  min_S_eigenvalue_hole_ = es_h.eigenvalues()(0);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_e(S_elec);
  min_S_eigenvalue_elec_ = es_e.eigenvalues()(0);

  XTP_LOG(Log::info, *pLog_)
      << "Smallest S eigenvalue hole/elec: " << min_S_eigenvalue_hole_ << " / "
      << min_S_eigenvalue_elec_ << flush;

  // -----------------------------------------------------------------------
  // Monomer MO energies — KS and QPpert (if available)
  //
  // KS energies: always present, directly from mos_.eigenvalues().
  //
  // QPpert energies: present only when a GW calculation was run on the
  // monomer. QPpert_energies_ is indexed over the GW window [qpmin, qpmax],
  // which must contain the DFT coupling range [Range_orbX.first, ..] for
  // the QP energies to be valid for our MOs. We verify this and fall back
  // to empty vectors (written as absent from the XML) if the window does
  // not cover the required range.
  // -----------------------------------------------------------------------

  // Helper lambda: extract QPpert energies for a given MO range,
  // returning an empty vector if GW is unavailable or the window
  // does not cover the requested range.
  auto ExtractQPEnergies = [&](const Orbitals& orb, Index mo_start,
                               Index n_mos) -> Eigen::VectorXd {
    if (!orb.hasQPpert()) return Eigen::VectorXd{};
    Index qpmin = orb.getGWAmin();
    Index qpmax = orb.getGWAmax();
    if (mo_start < qpmin || mo_start + n_mos - 1 > qpmax) {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp() << " WARNING: QPpert window [" << qpmin << "," << qpmax
          << "] does not fully cover MO range [" << mo_start << ","
          << mo_start + n_mos - 1 << "] -- QP energies omitted for this range"
          << flush;
      return Eigen::VectorXd{};
    }
    // QPpert_energies_ is indexed from qpmin, so offset accordingly
    return orb.QPpertEnergies().segment(mo_start - qpmin, n_mos) *
           votca::tools::conv::hrt2ev;
  };

  // KS energies — hole MOs: Range_orbA.first .. homoA
  moEnergiesA_hole_KS_ =
      orbitalsA.MOs().eigenvalues().segment(Range_orbA.first, holesA) *
      votca::tools::conv::hrt2ev;
  moEnergiesB_hole_KS_ =
      orbitalsB.MOs().eigenvalues().segment(Range_orbB.first, holesB) *
      votca::tools::conv::hrt2ev;

  // KS energies — electron MOs: lumo .. lumo + elecsX - 1
  moEnergiesA_elec_KS_ =
      orbitalsA.MOs().eigenvalues().segment(orbitalsA.getLumo(), elecsA) *
      votca::tools::conv::hrt2ev;
  moEnergiesB_elec_KS_ =
      orbitalsB.MOs().eigenvalues().segment(orbitalsB.getLumo(), elecsB) *
      votca::tools::conv::hrt2ev;

  // QPpert energies — extracted over the same ranges; empty if unavailable
  moEnergiesA_hole_QP_ = ExtractQPEnergies(orbitalsA, Range_orbA.first, holesA);
  moEnergiesB_hole_QP_ = ExtractQPEnergies(orbitalsB, Range_orbB.first, holesB);
  moEnergiesA_elec_QP_ =
      ExtractQPEnergies(orbitalsA, orbitalsA.getLumo(), elecsA);
  moEnergiesB_elec_QP_ =
      ExtractQPEnergies(orbitalsB, orbitalsB.getLumo(), elecsB);

  if (orbitalsA.hasQPpert()) {
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " QPpert energies available for fragment A" << flush;
  }
  if (orbitalsB.hasQPpert()) {
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " QPpert energies available for fragment B" << flush;
  }

  XTP_LOG(Log::error, *pLog_) << "Done with electronic couplings" << flush;
}

}  // namespace xtp
}  // namespace votca
