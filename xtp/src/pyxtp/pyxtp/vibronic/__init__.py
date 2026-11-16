"""Vibronic line shapes from normal modes and energy-gap gradients.

Workflow:

1. Hessian of the initial state at its minimum:
   ``FiniteDifferenceHessian`` (analytic forces, e.g. pyxtp DFT via
   ``xtp_ground_state_factory``) or ``read_orca_hessian``.
2. ``NormalModes.from_hessian``: frequencies and modes, translations and
   rotations projected out.
3. Gradient of the energy gap at the same geometry (BSE excitation energy,
   QP level, ...): ``couplings_from_gradient`` / ``couplings_from_forces``
   give lambda_nu and Huang-Rhys factors S_nu per mode.
4. ``vibronic_lineshape``: thermal Franck-Condon line shape (FFT of the
   generating function, quantum modes + classical part + broadening).
5. ``ensemble_spectrum``: sum over the vertical energies of a morphology.
6. Several levels at once (UPS, IPES): ``QPLevels.from_orb`` (dftgwbse
   option ``qp_gradients``) and ``level_spectrum``; figures with
   ``pyxtp.vibronic.plot`` (needs matplotlib).
"""

from .coupling import (VibronicCouplings, couplings_from_forces,
                       couplings_from_geometry, couplings_from_gradient,
                       read_gap_gradient)
from .hessian import (Displacement, FiniteDifferenceHessian,
                      xtp_ground_state_factory)
from .lineshape import LineShape, ensemble_spectrum, vibronic_lineshape
from .spectrum import QPLevels, Spectrum, add_gradients, level_spectrum
from .modes import (EV2INVCM, HBAR_OMEGA_EV, NormalModes, read_orca_engrad,
                    read_orca_hessian)

__all__ = ["VibronicCouplings", "couplings_from_forces",
           "couplings_from_geometry", "couplings_from_gradient",
           "read_gap_gradient",
           "Displacement", "FiniteDifferenceHessian",
           "xtp_ground_state_factory", "LineShape", "ensemble_spectrum",
           "vibronic_lineshape", "QPLevels", "Spectrum", "add_gradients",
           "level_spectrum",
           "EV2INVCM", "HBAR_OMEGA_EV", "NormalModes",
           "read_orca_engrad", "read_orca_hessian"]
