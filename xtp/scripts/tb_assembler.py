"""
TB assembler and solver for GW-BSE-DIPRO output.

Reads pairwise excitoncoupling XML files, assembles a sparse TB Hamiltonian
in the adaptive FE+CT basis (explicit CT states for pairs where downfolding
is unsafe, effective couplings otherwise), and computes eigenstates and
basic observables.

Site energy priority:
  1. External environment-corrected values passed in by caller (best).
  2. Monomer GW-BSE energies read from pair XML + optional electrostatic
     correction dict (pairwise-consistent, no environmental dressing).
  3. Averaged dimer diagonal H_FE_FE[i,i] (fallback with warning —
     pairwise-inconsistent, not recommended for production).

Usage:
    python tb_assembler.py
"""

import numpy as np
import xml.etree.ElementTree as ET
import scipy.linalg as la
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import warnings


# =============================================================================
# Data structures
# =============================================================================

@dataclass
class Diagnostics:
    xi: float
    pt_rm_discrepancy_eV: float
    downfolding_safe: bool


@dataclass
class PairData:
    """All data read from one excitoncoupling XML file."""
    frag_A: int
    frag_B: int
    n_FE: int           # total FE states in dimer (levA + levB)
    n_FE_A: int         # FE states on fragment A
    n_FE_B: int         # FE states on fragment B
    n_CT: int
    diagnostics: Diagnostics
    # monomer excitation energies [eV] — pairwise-consistent site energies
    monomer_energies_A: Optional[np.ndarray]   # shape (n_FE_A,)
    monomer_energies_B: Optional[np.ndarray]   # shape (n_FE_B,)
    # raw TB matrices (H in eV, S dimensionless)
    H_FE_FE: np.ndarray
    S_FE_FE: np.ndarray
    H_FE_CT: Optional[np.ndarray]
    S_FE_CT: Optional[np.ndarray]
    H_CT_CT: Optional[np.ndarray]
    S_CT_CT: Optional[np.ndarray]
    # effective couplings [eV]
    J_diag: np.ndarray   # shape (n_FE_A, n_FE_B)
    J_pert: np.ndarray   # shape (n_FE_A, n_FE_B)


# =============================================================================
# XML parsing
# =============================================================================

def parse_matrix(node) -> np.ndarray:
    """Parse a matrix node with row_N attributes into a numpy array."""
    rows = int(node.get("rows"))
    cols = int(node.get("cols"))
    mat = np.zeros((rows, cols))
    for i in range(rows):
        vals = [float(x) for x in node.get(f"row_{i}").split()]
        mat[i, :] = vals
    return mat


def parse_monomer_energies(mono_node, fragment: str) -> Optional[np.ndarray]:
    """
    Parse monomer_energies/fragmentA (or fragmentB) node.
    Returns array of energies in eV, or None if node absent.
    """
    frag_node = mono_node.find(fragment)
    if frag_node is None:
        return None
    energies = [float(e.get("eV")) for e in frag_node.findall("energy")]
    return np.array(energies)


def parse_coupling_xml(filepath: str, frag_A: int, frag_B: int) -> PairData:
    """
    Parse one excitoncoupling XML file.

    frag_A, frag_B: fragment indices (not stored in XML, supplied by caller).
    n_FE_A and n_FE_B are inferred from the coupling element state labels.
    """
    tree = ET.parse(filepath)
    root = tree.getroot()
    singlet = root.find(".//singlet")

    # --- effective couplings ---
    couplings = singlet.findall("coupling")
    stateA_labels = sorted(set(c.get("stateA") for c in couplings))
    stateB_labels = sorted(set(c.get("stateB") for c in couplings))
    n_FE_A = len(stateA_labels)
    n_FE_B = len(stateB_labels)
    stateA_idx = {s: i for i, s in enumerate(stateA_labels)}
    stateB_idx = {s: i for i, s in enumerate(stateB_labels)}

    J_diag = np.zeros((n_FE_A, n_FE_B))
    J_pert = np.zeros((n_FE_A, n_FE_B))
    for c in couplings:
        i = stateA_idx[c.get("stateA")]
        j = stateB_idx[c.get("stateB")]
        J_diag[i, j] = float(c.get("j_diag"))
        J_pert[i, j] = float(c.get("j_pert"))

    # --- monomer site energies ---
    mono_node = singlet.find("monomer_energies")
    if mono_node is not None:
        monomer_energies_A = parse_monomer_energies(mono_node, "fragmentA")
        monomer_energies_B = parse_monomer_energies(mono_node, "fragmentB")
    else:
        monomer_energies_A = None
        monomer_energies_B = None

    # --- diagnostics ---
    diag_node = singlet.find("diagnostics")
    diag = Diagnostics(
        xi=float(diag_node.get("xi")),
        pt_rm_discrepancy_eV=float(diag_node.get("pt_rm_discrepancy_eV")),
        downfolding_safe=(diag_node.get("downfolding_safe") == "true")
    )

    # --- TB matrices ---
    tb = singlet.find("tb_matrices")
    n_FE = int(tb.get("n_FE"))
    n_CT = int(tb.get("n_CT"))

    H_FE_FE = parse_matrix(tb.find("H_FE_FE"))
    S_FE_FE = parse_matrix(tb.find("S_FE_FE"))
    H_FE_CT = parse_matrix(tb.find("H_FE_CT")) if n_CT > 0 else None
    S_FE_CT = parse_matrix(tb.find("S_FE_CT")) if n_CT > 0 else None
    H_CT_CT = parse_matrix(tb.find("H_CT_CT")) if n_CT > 0 else None
    S_CT_CT = parse_matrix(tb.find("S_CT_CT")) if n_CT > 0 else None

    return PairData(
        frag_A=frag_A, frag_B=frag_B,
        n_FE=n_FE, n_FE_A=n_FE_A, n_FE_B=n_FE_B,
        n_CT=n_CT, diagnostics=diag,
        monomer_energies_A=monomer_energies_A,
        monomer_energies_B=monomer_energies_B,
        H_FE_FE=H_FE_FE, S_FE_FE=S_FE_FE,
        H_FE_CT=H_FE_CT, S_FE_CT=S_FE_CT,
        H_CT_CT=H_CT_CT, S_CT_CT=S_CT_CT,
        J_diag=J_diag, J_pert=J_pert
    )


# =============================================================================
# TB assembler
# =============================================================================

class TBAssembler:
    """
    Assembles the full TB Hamiltonian H and overlap S from pairwise dimer data.

    Basis ordering:
        [ FE_frag1_s1, FE_frag1_s2, ...,   <- FE block, ordered by fragment
          FE_frag2_s1, FE_frag2_s2, ...,
          ...
          CT_pair_AB_1, CT_pair_AB_2, ...,  <- CT block, explicit pairs only
          CT_pair_BC_1, ... ]

    Site energy precedence (FE diagonal):
        1. external_site_energies dict if provided
        2. monomer_energies from XML + optional electrostatic_corrections
        3. averaged dimer H_FE_FE diagonal (fallback, warns)

    CT diagonal always comes from H_CT_CT in the pair XML — no external
    correction needed since CT energies are intrinsically pairwise quantities.
    """

    def __init__(self, pairs: List[PairData], fragments: List[int]):
        self.pairs = pairs
        self.fragments = fragments
        self.FE_idx: Dict[Tuple[int, int], int] = {}
        self.CT_idx: Dict[Tuple[int, int, int], int] = {}
        self.explicit_CT_pairs: List[PairData] = []
        self.downfolded_pairs: List[PairData] = []
        self.total_dim: int = 0
        self._assign_indices()

    def _assign_indices(self):
        idx = 0
        n_FE_per_frag: Dict[int, int] = {}
        for pair in self.pairs:
            n_FE_per_frag[pair.frag_A] = pair.n_FE_A
            n_FE_per_frag[pair.frag_B] = pair.n_FE_B

        for frag in self.fragments:
            for i in range(n_FE_per_frag[frag]):
                self.FE_idx[(frag, i)] = idx
                idx += 1

        for pair in self.pairs:
            if not pair.diagnostics.downfolding_safe:
                self.explicit_CT_pairs.append(pair)
                for k in range(pair.n_CT):
                    self.CT_idx[(pair.frag_A, pair.frag_B, k)] = idx
                    idx += 1
            else:
                self.downfolded_pairs.append(pair)

        self.total_dim = idx

    def _resolve_site_energies(
            self,
            external_site_energies: Optional[Dict[Tuple[int,int], float]],
            electrostatic_corrections: Optional[Dict[Tuple[int,int], float]]
    ) -> Tuple[Dict[Tuple[int,int], float], str]:
        """
        Determine FE site energies using the priority hierarchy.
        Returns (site_energies dict, description of source used).
        """
        if external_site_energies is not None:
            for key in self.FE_idx:
                if key not in external_site_energies:
                    raise ValueError(
                        f"external_site_energies missing entry for "
                        f"fragment {key[0]} state {key[1]}")
            return dict(external_site_energies), "external (environment-corrected)"

        # Check whether monomer energies are present in all pair XMLs
        has_monomer_energies = all(
            pair.monomer_energies_A is not None and
            pair.monomer_energies_B is not None
            for pair in self.pairs
        )

        if has_monomer_energies:
            # Collect and verify consistency across pairs
            collected: Dict[Tuple[int,int], List[float]] = {}
            for pair in self.pairs:
                for i, eps in enumerate(pair.monomer_energies_A):
                    collected.setdefault((pair.frag_A, i), []).append(eps)
                for j, eps in enumerate(pair.monomer_energies_B):
                    collected.setdefault((pair.frag_B, j), []).append(eps)

            site_energies = {}
            for key, vals in collected.items():
                spread = max(vals) - min(vals)
                if spread > 1e-6:
                    warnings.warn(
                        f"Monomer energy for fragment {key[0]} state {key[1]} "
                        f"varies by {spread*1000:.3f} meV across pairs — "
                        f"check input files.", UserWarning)
                eps = float(np.mean(vals))
                if electrostatic_corrections is not None:
                    eps += electrostatic_corrections.get(key, 0.0)
                site_energies[key] = eps

            source = "monomer GW-BSE from XML"
            if electrostatic_corrections is not None:
                source += " + electrostatic corrections"
            return site_energies, source

        # Fallback: averaged dimer diagonal
        warnings.warn(
            "monomer_energies not found in XML. Falling back to averaged "
            "dimer diagonal H_FE_FE[i,i] as site energies. These are "
            "pairwise-inconsistent. Rerun with updated bsecoupling.cc to "
            "get monomer energies in the XML output.",
            UserWarning, stacklevel=3)

        energy_sum: Dict[Tuple[int,int], float] = {}
        energy_count: Dict[Tuple[int,int], int] = {}
        for pair in self.pairs:
            n_A = pair.n_FE_A
            for i in range(pair.n_FE_A):
                key = (pair.frag_A, i)
                energy_sum[key] = energy_sum.get(key, 0.0) + pair.H_FE_FE[i, i]
                energy_count[key] = energy_count.get(key, 0) + 1
            for j in range(pair.n_FE_B):
                key = (pair.frag_B, j)
                energy_sum[key] = energy_sum.get(key, 0.0) + \
                                  pair.H_FE_FE[n_A + j, n_A + j]
                energy_count[key] = energy_count.get(key, 0) + 1

        return ({k: energy_sum[k] / energy_count[k] for k in energy_sum},
                "averaged dimer diagonal (fallback, inconsistent)")

    def assemble(
            self,
            external_site_energies: Optional[Dict[Tuple[int,int], float]] = None,
            electrostatic_corrections: Optional[Dict[Tuple[int,int], float]] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Assemble and return (H, S) as dense matrices.

        Parameters
        ----------
        external_site_energies : dict {(frag_id, state_idx): energy_eV}
            Fully environment-corrected FE site energies (highest priority).
        electrostatic_corrections : dict {(frag_id, state_idx): delta_eV}
            Environmental shift to add on top of monomer XML energies.
            Only used when external_site_energies is None and monomer energies
            are present in the XML.
        """
        dim = self.total_dim
        H = np.zeros((dim, dim))
        S = np.zeros((dim, dim))
        np.fill_diagonal(S, 1.0)

        # --- FE site energies ---
        site_energies, source = self._resolve_site_energies(
            external_site_energies, electrostatic_corrections)
        self._site_energy_source = source  # store for reporting
        for (frag, i), idx in self.FE_idx.items():
            H[idx, idx] = site_energies[(frag, i)]

        # --- CT site energies and intra-pair CT-CT off-diagonal ---
        # Always from dimer H_CT_CT — intrinsically pairwise, no correction.
        for pair in self.explicit_CT_pairs:
            for k in range(pair.n_CT):
                kc = self.CT_idx[(pair.frag_A, pair.frag_B, k)]
                H[kc, kc] = pair.H_CT_CT[k, k]
                for l in range(pair.n_CT):
                    if l != k:
                        lc = self.CT_idx[(pair.frag_A, pair.frag_B, l)]
                        H[kc, lc] = pair.H_CT_CT[k, l]
                        S[kc, lc] = pair.S_CT_CT[k, l]

        # --- off-diagonal FE-FE blocks ---
        for pair in self.pairs:
            n_A = pair.n_FE_A
            idx_A = [self.FE_idx[(pair.frag_A, i)] for i in range(pair.n_FE_A)]
            idx_B = [self.FE_idx[(pair.frag_B, j)] for j in range(pair.n_FE_B)]

            if pair.diagnostics.downfolding_safe:
                # CT effects downfolded into J_diag
                for i, ia in enumerate(idx_A):
                    for j, jb in enumerate(idx_B):
                        H[ia, jb] = pair.J_diag[i, j]
                        H[jb, ia] = pair.J_diag[i, j]
                        S[ia, jb] = pair.S_FE_FE[i, n_A + j]
                        S[jb, ia] = pair.S_FE_FE[n_A + j, i]
            else:
                # Direct FE-FE coupling; CT states explicit in basis
                for i, ia in enumerate(idx_A):
                    for j, jb in enumerate(idx_B):
                        H[ia, jb] = pair.H_FE_FE[i, n_A + j]
                        H[jb, ia] = pair.H_FE_FE[n_A + j, i]
                        S[ia, jb] = pair.S_FE_FE[i, n_A + j]
                        S[jb, ia] = pair.S_FE_FE[n_A + j, i]
                # FE-CT couplings
                for i, ia in enumerate(idx_A):
                    for k in range(pair.n_CT):
                        kc = self.CT_idx[(pair.frag_A, pair.frag_B, k)]
                        H[ia, kc] = pair.H_FE_CT[i, k]
                        H[kc, ia] = pair.H_FE_CT[i, k]
                        S[ia, kc] = pair.S_FE_CT[i, k]
                        S[kc, ia] = pair.S_FE_CT[i, k]
                for j, jb in enumerate(idx_B):
                    for k in range(pair.n_CT):
                        kc = self.CT_idx[(pair.frag_A, pair.frag_B, k)]
                        H[jb, kc] = pair.H_FE_CT[n_A + j, k]
                        H[kc, jb] = pair.H_FE_CT[n_A + j, k]
                        S[jb, kc] = pair.S_FE_CT[n_A + j, k]
                        S[kc, jb] = pair.S_FE_CT[n_A + j, k]

        H = 0.5 * (H + H.T)
        S = 0.5 * (S + S.T)
        return H, S

    def basis_labels(self) -> List[str]:
        labels = [''] * self.total_dim
        for (frag, i), idx in self.FE_idx.items():
            labels[idx] = f"FE_f{frag}_s{i+1}"
        for (fA, fB, k), idx in self.CT_idx.items():
            labels[idx] = f"CT_f{fA}{fB}_k{k+1}"
        return labels


# =============================================================================
# TB solver
# =============================================================================

class TBSolver:
    """Solve the generalized eigenvalue problem H*c = E*S*c."""

    def __init__(self, H: np.ndarray, S: np.ndarray):
        self.H = H
        self.S = S
        self.evals: Optional[np.ndarray] = None
        self.evecs: Optional[np.ndarray] = None

    def solve(self) -> Tuple[np.ndarray, np.ndarray]:
        evals_S = np.linalg.eigvalsh(self.S)
        if evals_S.min() < 1e-6:
            warnings.warn(f"S has small eigenvalue {evals_S.min():.3e} — "
                          "basis may be near-linearly-dependent")
        self.evals, self.evecs = la.eigh(self.H, self.S)
        return self.evals, self.evecs

    def participation_ratio(self, evec: np.ndarray) -> float:
        p = np.abs(evec)**2
        return np.sum(p)**2 / np.sum(p**2)

    def fe_ct_weights(self, evec: np.ndarray,
                      assembler: TBAssembler) -> Tuple[float, float]:
        fe_w = sum(evec[idx]**2 for idx in assembler.FE_idx.values())
        ct_w = sum(evec[idx]**2 for idx in assembler.CT_idx.values())
        return fe_w, ct_w

    def fragment_weights(self, evec: np.ndarray,
                         assembler: TBAssembler,
                         fragments: List[int]) -> Dict[int, float]:
        return {
            frag: sum(evec[idx]**2
                      for (f, i), idx in assembler.FE_idx.items()
                      if f == frag)
            for frag in fragments
        }


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":

    fragments = [1, 2, 3]

    pair_files = [
        ("pair_1_2_excitoncoupling.xml", 1, 2),
        ("pair_1_3_excitoncoupling.xml", 1, 3),
        ("pair_2_3_excitoncoupling.xml", 2, 3),
    ]

    # --- parse ---
    print("=" * 60)
    print("Loading pairwise coupling data")
    print("=" * 60)
    pairs = []
    for filepath, fA, fB in pair_files:
        pair = parse_coupling_xml(filepath, fA, fB)
        pairs.append(pair)
        safe = "SAFE" if pair.diagnostics.downfolding_safe else "UNSAFE"
        print(f"\nPair {fA}-{fB}: downfolding {safe}")
        print(f"  xi = {pair.diagnostics.xi:.4f}  "
              f"PT/RM discrepancy = "
              f"{pair.diagnostics.pt_rm_discrepancy_eV*1000:.2f} meV")
        if pair.monomer_energies_A is not None:
            print(f"  Monomer energies A [eV]: {pair.monomer_energies_A}")
            print(f"  Monomer energies B [eV]: {pair.monomer_energies_B}")
        else:
            print(f"  Monomer energies: not in XML (old format)")
        if not pair.diagnostics.downfolding_safe:
            print(f"  -> CT states kept EXPLICIT")
            print(f"     CT energies [eV]: {np.diag(pair.H_CT_CT)}")
            print(f"     Max |H_FE_CT| = {np.max(np.abs(pair.H_FE_CT)):.4f} eV")
        else:
            print(f"  -> Using effective J_diag [eV]:")
            print(f"     {pair.J_diag}")

    # --- assemble ---
    print("\n" + "=" * 60)
    print("Assembling TB Hamiltonian")
    print("=" * 60)
    assembler = TBAssembler(pairs, fragments)
    labels = assembler.basis_labels()

    print(f"\nBasis: {assembler.total_dim} states total  "
          f"({len(assembler.FE_idx)} FE + {len(assembler.CT_idx)} CT)")
    for i, label in enumerate(labels):
        print(f"  [{i}] {label}")

    # Assemble — uses monomer XML energies as site energies.
    # To add electrostatic corrections per fragment/state:
    #   electrostatic_corrections = {(1, 0): +0.05, (1, 1): +0.03, ...}
    # To use fully external site energies:
    #   external_site_energies = {(1, 0): 3.80, ...}
    H, S = assembler.assemble()

    print(f"\nSite energy source: {assembler._site_energy_source}")
    print(f"\nFE site energies used [eV]:")
    for (frag, i), idx in sorted(assembler.FE_idx.items()):
        print(f"  Fragment {frag} state {i+1}: {H[idx, idx]:.6f}")
    print(f"\nCT site energies [eV]:")
    for (fA, fB, k), idx in sorted(assembler.CT_idx.items()):
        print(f"  CT f{fA}-f{fB} k{k+1}: {H[idx, idx]:.6f}")

    evals_S = np.linalg.eigvalsh(S)
    print(f"\nS eigenvalues: min={evals_S.min():.6f}  max={evals_S.max():.6f}")
    print(f"S positive definite: {np.all(evals_S > 0)}")

    # --- solve ---
    print("\n" + "=" * 60)
    print("TB eigenstates")
    print("=" * 60)
    solver = TBSolver(H, S)
    evals, evecs = solver.solve()

    print(f"\n{'#':>4}  {'E (eV)':>10}  {'IPR':>6}  "
          f"{'FE':>6}  {'CT':>6}  "
          f"{'f1':>6}  {'f2':>6}  {'f3':>6}  Dominant")
    print("-" * 75)
    for n in range(assembler.total_dim):
        evec = evecs[:, n]
        ipr = solver.participation_ratio(evec)
        fe_w, ct_w = solver.fe_ct_weights(evec, assembler)
        fw = solver.fragment_weights(evec, assembler, fragments)
        dom = labels[np.argmax(np.abs(evec)**2)]
        print(f"{n+1:>4}  {evals[n]:>10.4f}  {ipr:>6.2f}  "
              f"{fe_w:>6.3f}  {ct_w:>6.3f}  "
              f"{fw[1]:>6.3f}  {fw[2]:>6.3f}  {fw[3]:>6.3f}  {dom}")

    # --- compare site energy sources ---
    print("\n" + "=" * 60)
    print("Site energy source comparison (FE states only)")
    print("=" * 60)
    n_FE = len(assembler.FE_idx)

    # Fallback: force dimer diagonal by temporarily hiding monomer energies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        saved = [(p.monomer_energies_A, p.monomer_energies_B) for p in pairs]
        for p in pairs:
            p.monomer_energies_A = None
            p.monomer_energies_B = None
        H_fb, S_fb = assembler.assemble()
        for p, (a, b) in zip(pairs, saved):
            p.monomer_energies_A = a
            p.monomer_energies_B = b
    evals_fb, _ = la.eigh(H_fb, S_fb)

    print(f"\n{'#':>4}  {'E_monomer (eV)':>16}  "
          f"{'E_dimer_avg (eV)':>17}  {'Delta (meV)':>12}")
    print("-" * 55)
    for n in range(n_FE):
        dE = (evals[n] - evals_fb[n]) * 1000
        print(f"{n+1:>4}  {evals[n]:>16.4f}  "
              f"{evals_fb[n]:>17.4f}  {dE:>12.2f}")

    # --- eigenvector composition ---
    print("\n" + "=" * 60)
    print("Eigenvector composition (|coeff| > 0.1)")
    print("=" * 60)
    for n in range(assembler.total_dim):
        evec = evecs[:, n]
        print(f"\nState {n+1}  E = {evals[n]:.4f} eV:")
        for c, label in zip(evec, labels):
            if abs(c) > 0.1:
                print(f"  {c:+.4f}  {label}")
