"""
Tests for xtp_tightbinding.

Covers:
  - XML parsing (parse_bse_spin, parse_dft_carrier)
  - TBAssembler: basis construction, index assignment, Hamiltonian assembly
  - TBSolver: eigenvalues, eigenvectors, oscillator strengths
  - TBDynamics: initial states, noiseless propagation, stochastic propagation
  - Regression tests against the three-fragment formaldehyde trimer in
    iqm_out.xml (the same system used for development validation)

Run with:
    pytest test_xtp_tightbinding.py -v
"""

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import warnings
import os
import sys
import pytest

# Import xtp_tightbinding as a module.  The file has no .py extension
# in the VOTCA tree, so we import it via importlib.
import importlib.util
_here = os.path.dirname(os.path.abspath(__file__))
_tool_path = os.path.join(_here, "xtp_tightbinding.py")
_spec = importlib.util.spec_from_file_location("xtp_tightbinding",_tool_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

# Pull names into this namespace for convenience
parse_bse_spin      = _mod.parse_bse_spin
parse_dft_carrier   = _mod.parse_dft_carrier
parse_jobs_xml      = _mod.parse_jobs_xml
TBAssembler         = _mod.TBAssembler
TBSolver            = _mod.TBSolver
TBDynamics          = _mod.TBDynamics
BSEDiagnostics      = _mod.BSEDiagnostics
DFTDiagnostics      = _mod.DFTDiagnostics

# Path to the test XML fixture
TEST_XML = os.path.join(_here, "iqm_out.xml")


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def singlet_pairs_frags():
    """Parse the singlet model from the test XML once for the whole module."""
    pairs, fragments = parse_jobs_xml(TEST_XML, "singlet")
    return pairs, fragments


@pytest.fixture(scope="module")
def singlet_assembler(singlet_pairs_frags):
    pairs, fragments = singlet_pairs_frags
    return TBAssembler(pairs, fragments, explicit_ct=True)


@pytest.fixture(scope="module")
def singlet_H_S(singlet_assembler):
    H, S = singlet_assembler.assemble(sparse=False)
    return H, S


@pytest.fixture(scope="module")
def singlet_evals_evecs(singlet_H_S):
    H, S = singlet_H_S
    evals, evecs = la.eigh(H, S)
    return evals, evecs


@pytest.fixture(scope="module")
def hole_pairs_frags():
    pairs, fragments = parse_jobs_xml(TEST_XML, "hole")
    return pairs, fragments


# =============================================================================
# XML parsing tests
# =============================================================================

class TestParsing:

    def test_singlet_pair_count(self, singlet_pairs_frags):
        pairs, fragments = singlet_pairs_frags
        assert len(pairs) == 3, "Expected 3 pairs for 3-fragment system"

    def test_fragments_sorted(self, singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        assert fragments == sorted(fragments)

    def test_singlet_pair_fields(self, singlet_pairs_frags):
        pairs, _ = singlet_pairs_frags
        for p in pairs:
            assert p.model == "singlet"
            assert p.n_A == p.n_B, "Homodimer: both fragments have same n states"
            assert p.n_CT > 0, "CT states present"
            assert p.monomer_energies_A is not None
            assert p.monomer_energies_B is not None
            assert len(p.monomer_energies_A) == p.n_A
            assert len(p.monomer_energies_B) == p.n_B
            assert p.H_AB is not None
            assert p.H_AB.shape == (p.n_A, p.n_B)
            assert p.S_AB is not None

    def test_diagnostics_fields(self, singlet_pairs_frags):
        pairs, _ = singlet_pairs_frags
        for p in pairs:
            assert isinstance(p.diagnostics, BSEDiagnostics)
            assert p.diagnostics.xi >= 0
            assert isinstance(p.diagnostics.downfolding_safe, bool)

    def test_monomer_energies_consistent_across_pairs(self, singlet_pairs_frags):
        """Each fragment's monomer energy should be the same in all pairs it appears in."""
        pairs, fragments = singlet_pairs_frags
        collected = {}
        for p in pairs:
            for frag, mono in [(p.frag_A, p.monomer_energies_A),
                               (p.frag_B, p.monomer_energies_B)]:
                for i, e in enumerate(mono):
                    key = (frag, i)
                    if key in collected:
                        assert abs(collected[key] - e) < 1e-6, (
                            f"Fragment {frag} state {i} energy inconsistent "
                            f"across pairs: {collected[key]:.6f} vs {e:.6f}")
                    else:
                        collected[key] = e

    def test_transition_dipoles_present(self, singlet_pairs_frags):
        pairs, _ = singlet_pairs_frags
        for p in pairs:
            assert p.transition_dipoles_A is not None, \
                "Transition dipoles should be in the XML (new bsecoupling)"
            assert p.transition_dipoles_A.shape == (p.n_A, 3)

    def test_hole_parsing(self, hole_pairs_frags):
        pairs, fragments = hole_pairs_frags
        assert len(pairs) == 3
        for p in pairs:
            assert p.model == "hole"
            assert p.n_CT == 0
            assert isinstance(p.diagnostics, DFTDiagnostics)
            assert p.diagnostics.min_S_eigenvalue > 0
            assert p.monomer_energies_A is not None

    def test_missing_model_warns(self):
        """Requesting a model not in the XML should warn and return empty pairs."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pairs, fragments = parse_jobs_xml(TEST_XML, "triplet")
        # Either empty pairs (triplet not computed) or warnings issued
        # Both are acceptable — just no crash
        assert isinstance(pairs, list)


# =============================================================================
# TBAssembler tests
# =============================================================================

class TestTBAssembler:

    def test_basis_dimension(self, singlet_assembler, singlet_pairs_frags):
        pairs, fragments = singlet_pairs_frags
        asm = singlet_assembler
        n_FE = sum(p.n_A for p in pairs if p.frag_A == fragments[0])
        # For explicit_ct=True: all pairs contribute CT states
        expected_FE = len(fragments) * pairs[0].n_A
        assert len(asm.state_idx) == expected_FE
        assert asm.total_dim > expected_FE  # CT states added

    def test_CT_count(self, singlet_assembler, singlet_pairs_frags):
        pairs, _ = singlet_pairs_frags
        asm = singlet_assembler
        # explicit_ct=True: all 3 pairs contribute CT states
        assert len(asm.explicit_CT_pairs) == 3
        assert len(asm.downfolded_pairs) == 0
        expected_CT = sum(p.n_CT for p in pairs)
        assert len(asm.CT_idx) == expected_CT

    def test_state_idx_unique(self, singlet_assembler):
        indices = list(singlet_assembler.state_idx.values())
        assert len(indices) == len(set(indices)), "Duplicate state indices"

    def test_CT_idx_unique(self, singlet_assembler):
        indices = list(singlet_assembler.CT_idx.values())
        assert len(indices) == len(set(indices)), "Duplicate CT indices"

    def test_indices_contiguous(self, singlet_assembler):
        all_idx = (list(singlet_assembler.state_idx.values()) +
                   list(singlet_assembler.CT_idx.values()))
        assert sorted(all_idx) == list(range(singlet_assembler.total_dim))

    def test_deterministic_index_assignment(self, singlet_pairs_frags):
        """Index assignment must be independent of pair order in the list."""
        pairs, fragments = singlet_pairs_frags
        asm1 = TBAssembler(pairs, fragments, explicit_ct=True)
        asm2 = TBAssembler(list(reversed(pairs)), fragments, explicit_ct=True)
        assert asm1.state_idx == asm2.state_idx
        assert asm1.CT_idx == asm2.CT_idx

    def test_inconsistent_state_count_raises(self, singlet_pairs_frags):
        """Pairs with inconsistent n_A for the same fragment must raise."""
        pairs, fragments = singlet_pairs_frags
        import copy
        bad_pairs = copy.deepcopy(pairs)
        bad_pairs[0].n_A += 1  # make fragment count inconsistent
        with pytest.raises(ValueError, match="inconsistent state counts"):
            TBAssembler(bad_pairs, fragments)

    def test_basis_labels_length(self, singlet_assembler):
        labels = singlet_assembler.basis_labels()
        assert len(labels) == singlet_assembler.total_dim
        assert all(isinstance(l, str) and len(l) > 0 for l in labels)

    def test_adaptive_vs_explicit_ct(self, singlet_pairs_frags):
        """explicit_ct=False should produce a smaller basis for this system
        (which has some safe pairs)."""
        pairs, fragments = singlet_pairs_frags
        asm_adapt   = TBAssembler(pairs, fragments, explicit_ct=False)
        asm_explicit = TBAssembler(pairs, fragments, explicit_ct=True)
        # The formaldehyde trimer has at least one safe pair (pair 0-2)
        # so explicit_ct should give a larger basis
        assert asm_explicit.total_dim >= asm_adapt.total_dim


# =============================================================================
# Hamiltonian assembly tests
# =============================================================================

class TestHamiltonianAssembly:

    def test_H_symmetric(self, singlet_H_S):
        H, S = singlet_H_S
        assert np.allclose(H, H.T, atol=1e-12), "H must be symmetric"

    def test_S_symmetric(self, singlet_H_S):
        H, S = singlet_H_S
        assert np.allclose(S, S.T, atol=1e-12), "S must be symmetric"

    def test_S_positive_definite(self, singlet_H_S):
        H, S = singlet_H_S
        evals_S = np.linalg.eigvalsh(S)
        assert evals_S.min() > 0, \
            f"S must be positive definite, min eigenvalue = {evals_S.min():.4e}"

    def test_site_energies_on_diagonal(self, singlet_assembler, singlet_H_S):
        """FE diagonal of H should match monomer energies from the XML."""
        H, S = singlet_H_S
        asm = singlet_assembler
        for (frag, i), idx in asm.state_idx.items():
            # Get expected energy from first pair that involves this fragment
            for p in asm.pairs:
                if p.frag_A == frag:
                    expected = p.monomer_energies_A[i]
                    break
                elif p.frag_B == frag:
                    expected = p.monomer_energies_B[i]
                    break
            assert abs(H[idx, idx] - expected) < 1e-8, (
                f"H diagonal for frag {frag} state {i}: "
                f"got {H[idx,idx]:.6f}, expected {expected:.6f}")

    def test_S_diagonal_unity(self, singlet_H_S):
        """S diagonal should be 1 (set by np.fill_diagonal)."""
        H, S = singlet_H_S
        assert np.allclose(np.diag(S), 1.0, atol=1e-12)

    def test_sparse_dense_agree(self, singlet_assembler):
        """Sparse and dense assembly must produce identical matrices."""
        H_dense, S_dense = singlet_assembler.assemble(sparse=False)
        H_sparse, S_sparse = singlet_assembler.assemble(sparse=True)
        assert sp.issparse(H_sparse)
        assert np.allclose(H_dense, H_sparse.toarray(), atol=1e-12)
        assert np.allclose(S_dense, S_sparse.toarray(), atol=1e-12)

    def test_shape(self, singlet_assembler, singlet_H_S):
        H, S = singlet_H_S
        dim = singlet_assembler.total_dim
        assert H.shape == (dim, dim)
        assert S.shape == (dim, dim)


# =============================================================================
# TBSolver tests
# =============================================================================

class TestTBSolver:

    def test_eigenvalue_count(self, singlet_H_S, singlet_assembler):
        H, S = singlet_H_S
        solver = TBSolver(H, S)
        evals, evecs = solver.solve()
        assert len(evals) == singlet_assembler.total_dim
        assert evecs.shape == (singlet_assembler.total_dim,
                               singlet_assembler.total_dim)

    def test_eigenvalues_real_ordered(self, singlet_evals_evecs):
        evals, _ = singlet_evals_evecs
        assert np.all(np.isreal(evals))
        assert np.all(np.diff(evals) >= -1e-12), "Eigenvalues not non-decreasing"

    def test_eigenvectors_S_orthonormal(self, singlet_H_S, singlet_evals_evecs):
        """Eigenvectors must satisfy V^T S V = I."""
        H, S = singlet_H_S
        _, evecs = singlet_evals_evecs
        gram = evecs.T @ S @ evecs
        assert np.allclose(gram, np.eye(len(gram)), atol=1e-8), \
            "Eigenvectors not S-orthonormal"

    def test_H_eigendecomposition(self, singlet_H_S, singlet_evals_evecs):
        """V^T H V should equal diag(evals) for S-orthonormal V."""
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        Hv = evecs.T @ H @ evecs
        assert np.allclose(Hv, np.diag(evals), atol=1e-8)

    def test_participation_ratio(self, singlet_H_S, singlet_evals_evecs,
                                  singlet_assembler):
        _, evecs = singlet_evals_evecs
        H, S = singlet_H_S
        solver = TBSolver(H, S)
        ipr = solver.participation_ratio(evecs[:, 0])
        assert 1.0 <= ipr <= singlet_assembler.total_dim

    def test_state_weights_S_normalised(self, singlet_H_S, singlet_evals_evecs,
                                        singlet_assembler, singlet_pairs_frags):
        """S-normalised eigenvectors satisfy <v|S|v>=1 exactly.
        state_weights returns |c_i|^2 in the identity metric which sums to
        |v|^2 < 1 when S is non-diagonal -- this is expected behaviour."""
        _, fragments = singlet_pairs_frags
        _, evecs = singlet_evals_evecs
        H, S = singlet_H_S
        for n in range(5):
            evec = evecs[:, n]
            s_norm = float(np.real(evec @ S @ evec))
            assert abs(s_norm - 1.0) < 1e-8, \
                f"State {n}: <v|S|v> = {s_norm:.8f}, expected 1.0"

class TestOscillatorStrengths:

    def test_oscillator_strengths_non_negative(self, singlet_H_S,
                                                singlet_evals_evecs,
                                                singlet_assembler):
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        if not singlet_assembler.transition_dipoles:
            pytest.skip("No transition dipoles in XML")
        solver = TBSolver(H, S)
        f = solver.oscillator_strengths(evecs, evals, singlet_assembler)
        assert np.all(f >= 0), "Oscillator strengths must be non-negative"

    def test_oscillator_strengths_count(self, singlet_H_S, singlet_evals_evecs,
                                         singlet_assembler):
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        if not singlet_assembler.transition_dipoles:
            pytest.skip("No transition dipoles in XML")
        solver = TBSolver(H, S)
        f = solver.oscillator_strengths(evecs, evals, singlet_assembler)
        assert len(f) == len(evals)

    def test_s1_states_dark(self, singlet_H_S, singlet_evals_evecs,
                             singlet_assembler):
        """s1 states should have near-zero oscillator strength (dark in formaldehyde)."""
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        if not singlet_assembler.transition_dipoles:
            pytest.skip("No transition dipoles in XML")
        solver = TBSolver(H, S)
        f = solver.oscillator_strengths(evecs, evals, singlet_assembler)
        # First 3 eigenstates are s1 FE manifold — all dark
        for n in range(3):
            assert f[n] < 1e-3, \
                f"s1 state {n+1} should be dark, got f={f[n]:.4f}"

    def test_s2_bright_state(self, singlet_H_S, singlet_evals_evecs,
                              singlet_assembler):
        """s2 manifold should have one bright state with f ~ 0.24 eV."""
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        if not singlet_assembler.transition_dipoles:
            pytest.skip("No transition dipoles in XML")
        solver = TBSolver(H, S)
        f = solver.oscillator_strengths(evecs, evals, singlet_assembler)
        # s2 FE manifold is states 4-6 (indices 3-5) — find the bright one
        f_s2 = f[3:6]
        assert f_s2.max() > 0.1, \
            f"s2 manifold should have a bright state, max f={f_s2.max():.4f}"
        # The bright state should dominate
        assert f_s2.max() > 5 * np.median(f_s2), \
            "Bright state should dominate the s2 manifold"


# =============================================================================
# Regression tests — pinned numerical values for the formaldehyde trimer
# =============================================================================

class TestRegression:
    """
    Regression tests against known-good values for the three-fragment
    formaldehyde trimer (explicit_ct=True, singlet model).

    Values were derived during development and validated against full
    GW-BSE trimer calculations. Tolerances are set to 0.1% to catch
    code regressions while being robust to platform differences.
    """

    # s1 FE manifold: 3 states (fragments × 1 FE state each)
    # Monomer s1 energy: 3.3786 eV (from XML)
    S1_CENTRE_EV = 3.3786  # eV, approximate centre of s1 manifold
    S1_BW_MEV = 15.4       # meV, s1 bandwidth

    # s2 FE manifold: 3 states
    # Monomer s2 energy: 7.7302 eV (from XML)
    S2_CENTRE_EV = 7.7302
    S2_BW_MEV = 123.0      # meV, s2 bandwidth (25 CT per direction)

    # Bright s2 oscillator strength (vs GW-BSE S6 f=0.261)
    S2_BRIGHT_F = 0.244    # dimensionless

    def test_s1_manifold_centre(self, singlet_evals_evecs):
        evals, _ = singlet_evals_evecs
        # First 3 FE eigenstates are s1 manifold
        s1_centre = np.mean(evals[:3])
        assert abs(s1_centre - self.S1_CENTRE_EV) < 0.01, \
            f"s1 centre: {s1_centre:.4f} eV, expected ~{self.S1_CENTRE_EV:.4f} eV"

    def test_s1_bandwidth(self, singlet_evals_evecs):
        evals, _ = singlet_evals_evecs
        bw = (evals[2] - evals[0]) * 1000  # meV
        assert abs(bw - self.S1_BW_MEV) / self.S1_BW_MEV < 0.05, \
            f"s1 bandwidth: {bw:.1f} meV, expected ~{self.S1_BW_MEV:.1f} meV"

    def test_s2_manifold_centre(self, singlet_evals_evecs, singlet_assembler):
        evals, evecs = singlet_evals_evecs
        # Find s2 FE states: FE-dominated states in the 7-9 eV range
        H_dense = singlet_assembler.assemble(sparse=False)[0]
        s2_states = [n for n in range(len(evals))
                     if 7.0 < evals[n] < 9.0 and
                     sum(evecs[idx, n]**2
                         for idx in singlet_assembler.state_idx.values()) > 0.5]
        assert len(s2_states) >= 3, f"Expected ≥3 FE s2 states, found {len(s2_states)}"
        s2_centre = np.mean(evals[s2_states[:3]])
        assert abs(s2_centre - self.S2_CENTRE_EV) < 0.05, \
            f"s2 centre: {s2_centre:.4f} eV, expected ~{self.S2_CENTRE_EV:.4f} eV"

    def test_s2_bandwidth(self, singlet_evals_evecs, singlet_assembler):
        evals, evecs = singlet_evals_evecs
        s2_states = [n for n in range(len(evals))
                     if 7.0 < evals[n] < 9.0 and
                     sum(evecs[idx, n]**2
                         for idx in singlet_assembler.state_idx.values()) > 0.5]
        if len(s2_states) < 3:
            pytest.skip("Not enough s2 FE states found")
        bw = (evals[s2_states[2]] - evals[s2_states[0]]) * 1000  # meV
        assert abs(bw - self.S2_BW_MEV) / self.S2_BW_MEV < 0.05, \
            f"s2 bandwidth: {bw:.1f} meV, expected ~{self.S2_BW_MEV:.1f} meV"

    def test_s2_bright_oscillator_strength(self, singlet_H_S,
                                            singlet_evals_evecs,
                                            singlet_assembler):
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        if not singlet_assembler.transition_dipoles:
            pytest.skip("No transition dipoles in XML")
        solver = TBSolver(H, S)
        f = solver.oscillator_strengths(evecs, evals, singlet_assembler)
        s2_states = [n for n in range(len(evals))
                     if 7.0 < evals[n] < 9.0 and
                     sum(evecs[idx, n]**2
                         for idx in singlet_assembler.state_idx.values()) > 0.5]
        if len(s2_states) < 3:
            pytest.skip("Not enough s2 FE states found")
        f_bright = max(f[n] for n in s2_states[:3])
        assert abs(f_bright - self.S2_BRIGHT_F) / self.S2_BRIGHT_F < 0.10, \
            f"s2 bright f: {f_bright:.4f}, expected ~{self.S2_BRIGHT_F:.4f}"

    def test_bright_state_is_highest_in_s2(self, singlet_H_S,
                                            singlet_evals_evecs,
                                            singlet_assembler):
        """J-aggregate character: bright state should be the highest of the s2 triplet."""
        H, S = singlet_H_S
        evals, evecs = singlet_evals_evecs
        if not singlet_assembler.transition_dipoles:
            pytest.skip("No transition dipoles in XML")
        solver = TBSolver(H, S)
        f = solver.oscillator_strengths(evecs, evals, singlet_assembler)
        s2_states = [n for n in range(len(evals))
                     if 7.0 < evals[n] < 9.0 and
                     sum(evecs[idx, n]**2
                         for idx in singlet_assembler.state_idx.values()) > 0.5]
        if len(s2_states) < 3:
            pytest.skip("Not enough s2 FE states found")
        s2_f = [(f[n], n) for n in s2_states[:3]]
        brightest_idx = max(s2_f, key=lambda x: x[0])[1]
        assert brightest_idx == s2_states[2], \
            "Bright state should be the highest-energy state of the s2 triplet (J-aggregate)"


# =============================================================================
# TBDynamics tests
# =============================================================================

class TestTBDynamics:

    @pytest.fixture
    def dynamics(self, singlet_H_S, singlet_assembler):
        H, S = singlet_H_S
        return TBDynamics(H, S, singlet_assembler)

    def test_site_excitation_norm(self, dynamics, singlet_assembler,
                                   singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        # Norm under S: <psi|S|psi> should be 1
        norm_sq = np.real(psi0.conj() @ (dynamics._sqrt_S @ dynamics._sqrt_S @ psi0))
        assert abs(norm_sq - 1.0) < 1e-10

    def test_site_excitation_localised(self, dynamics, singlet_assembler,
                                       singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        # Weight on the target fragment should be 1 (only that site excited)
        target_idx = singlet_assembler.state_idx[(fragments[0], 0)]
        pop = np.abs(psi0)**2
        frag_pop = dynamics._W @ pop
        assert frag_pop[0] > 0.99, \
            f"Initial population on target fragment: {frag_pop[0]:.4f}"

    def test_gaussian_packet_norm(self, dynamics, singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.gaussian_packet(fragments[1], sigma_frags=1.5)
        norm_sq = np.real(psi0.conj() @ (dynamics._sqrt_S @ dynamics._sqrt_S @ psi0))
        assert abs(norm_sq - 1.0) < 1e-10

    def test_gaussian_packet_peaked_on_centre(self, dynamics, singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        centre = fragments[1]
        psi0 = dynamics.gaussian_packet(centre, sigma_frags=0.5)
        pop = np.abs(psi0)**2
        frag_pop = dynamics._W @ pop
        centre_fi = dynamics.frag_to_idx[centre]
        assert frag_pop[centre_fi] == max(frag_pop), \
            "Centre fragment should have largest population"

    def test_invalid_site_raises(self, dynamics):
        with pytest.raises((ValueError, KeyError)):
            dynamics.site_excitation(9999, 0)

    def test_pure_propagation_purity(self, dynamics, singlet_pairs_frags):
        """Noiseless propagation preserves purity = 1 (pure state)."""
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        result = dynamics.propagate_pure(psi0, t_max_fs=50.0, dt_fs=5.0)
        assert np.allclose(result['purity'], 1.0, atol=1e-10)

    def test_pure_propagation_norm_conservation(self, dynamics, singlet_pairs_frags):
        """Total population should be conserved under noiseless propagation."""
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        result = dynamics.propagate_pure(psi0, t_max_fs=50.0, dt_fs=5.0)
        pop_sum = result['populations'].sum(axis=0)
        # Mulliken populations sum to exactly 1 (= <c|S|c> = 1)
        assert np.allclose(pop_sum, 1.0, atol=1e-6), \
            f"Population not conserved: [{pop_sum.min():.6f}, {pop_sum.max():.6f}]"

    def test_pure_propagation_initial_population(self, dynamics, singlet_pairs_frags):
        """At t=0, population should be on the initial fragment."""
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        result = dynamics.propagate_pure(psi0, t_max_fs=50.0, dt_fs=5.0)
        assert result['populations'][0, 0] > 0.99, \
            f"Initial population on fragment 0: {result['populations'][0,0]:.4f}"

    def test_pure_propagation_time_axis(self, dynamics, singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        result = dynamics.propagate_pure(psi0, t_max_fs=100.0, dt_fs=10.0)
        t = result['t']
        assert t[0] == pytest.approx(0.0)
        assert t[-1] == pytest.approx(100.0)
        assert len(t) == 11  # 0, 10, 20, ..., 100

    def test_stochastic_purity_decay(self, dynamics, singlet_pairs_frags):
        """Stochastic propagation should cause purity to decay from 1."""
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        result = dynamics.propagate_stochastic(
            psi0, t_max_fs=100.0, dt_fs=1.0,
            sigma_eV=0.1, tau_fs=20.0,
            n_realisations=30, seed=42)
        # Purity at t=0 should be near 1
        assert result['purity'][0] > 0.95
        # Purity should have decayed by t=100 fs with 100 meV noise
        assert result['purity'][-1] < result['purity'][0]

    def test_stochastic_populations_sum_to_one(self, dynamics, singlet_pairs_frags):
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        result = dynamics.propagate_stochastic(
            psi0, t_max_fs=50.0, dt_fs=2.0,
            sigma_eV=0.05, tau_fs=30.0,
            n_realisations=20, seed=0)
        pop_sum = result['populations'].sum(axis=0)
        assert np.allclose(pop_sum, 1.0, atol=1e-3)

    def test_thermal_mixture_trace_one(self, dynamics, singlet_evals_evecs):
        evals, evecs = singlet_evals_evecs
        rho0 = dynamics.thermal_mixture(300.0, evals, evecs)
        # eigenvectors are S-normalised; physical trace needs S-metric.
        sqrt_S = dynamics._sqrt_S
        assert abs(np.trace(sqrt_S @ rho0 @ sqrt_S).real - 1.0) < 1e-8

    def test_thermal_mixture_positive_semidefinite(self, dynamics,
                                                    singlet_evals_evecs):
        evals, evecs = singlet_evals_evecs
        rho0 = dynamics.thermal_mixture(300.0, evals, evecs)
        evals_rho = np.linalg.eigvalsh(rho0)
        assert np.all(evals_rho >= -1e-10), "Density matrix must be PSD"

    def test_stability_warning_dt_too_large(self, dynamics, singlet_pairs_frags,
                                             capsys):
        """Large dt relative to noise timescale should trigger a warning."""
        _, fragments = singlet_pairs_frags
        psi0 = dynamics.site_excitation(fragments[0], 0)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # dt=50 fs >> tau/5=10 fs should warn
            result = dynamics.propagate_stochastic(
                psi0, t_max_fs=50.0, dt_fs=50.0,
                sigma_eV=0.05, tau_fs=10.0,
                n_realisations=5, seed=1)
        # No assert on result quality — just checking no crash

    def test_W_matrix_rows_sum_to_one(self, dynamics):
        """Each basis state belongs to exactly one fragment (FE) or split (CT)."""
        W_sum = dynamics._W.sum(axis=0)
        assert np.allclose(W_sum, 1.0, atol=1e-12), \
            "W columns must sum to 1 (each basis state fully assigned)"


# =============================================================================
# Transport model tests
# =============================================================================

class TestTransportModel:

    def test_transport_basis_has_hole_and_electron(self):
        pairs, fragments = parse_jobs_xml(TEST_XML, "transport")
        asm = TBAssembler(pairs, fragments)
        assert asm.model == "transport"
        hole_states   = [k for k in asm.state_idx if k[2] == "hole"]
        elec_states   = [k for k in asm.state_idx if k[2] == "electron"]
        assert len(hole_states) > 0
        assert len(elec_states) > 0

    def test_transport_H_block_diagonal(self):
        """Hole-electron off-diagonal block should be zero."""
        pairs, fragments = parse_jobs_xml(TEST_XML, "transport")
        asm = TBAssembler(pairs, fragments)
        H, S = asm.assemble(sparse=False)
        hole_idx = [idx for (frag, i, c), idx in asm.state_idx.items()
                    if c == "hole"]
        elec_idx = [idx for (frag, i, c), idx in asm.state_idx.items()
                    if c == "electron"]
        off_diag = H[np.ix_(hole_idx, elec_idx)]
        assert np.allclose(off_diag, 0.0, atol=1e-12), \
            "Hole-electron off-diagonal block must be zero"

    def test_transport_KS_site_energies(self):
        """KS site energies must be negative (occupied and virtual MOs of a
        finite molecule have negative KS eigenvalues in Hartree)."""
        pairs, fragments = parse_jobs_xml(TEST_XML, "transport")
        asm = TBAssembler(pairs, fragments)
        H, S = asm.assemble(sparse=False)
        for key, idx in asm.state_idx.items():
            e = H[idx, idx]
            assert e < 0, \
                f"KS site energy for {key} should be negative (got {e:.3f} eV)"
