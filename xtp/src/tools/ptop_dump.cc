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

// Dev/debug tool: NOT part of the production xtp_* tool set (no --help,
// manpage, or sphinx-doc infrastructure). Dumps the background shell
// (PolarTop::BGN()) of a legacy .ptop checkpoint (see
// PolarTop::SaveToDrive/LoadFromDrive, boost::serialization-based) as a
// CSV, one row per site, with positions, permanent fields, and induced
// dipoles converted into the same units the new EwaldRegistry/PolarSite
// code uses (bohr length, atomic-unit-style fields/dipoles), for direct
// comparison against hdf5_dump's own output on the same schema.
//
// Deliberately links ONLY against legacy xtp/ewald/ code (not the new
// EwaldRegistry machinery at all) to avoid any risk of symbol/namespace
// friction between APolarSite and PolarSite -- hdf5_dump is the
// equivalent tool for the new side, kept fully separate for the same
// reason; the two outputs are compared by a separate script, not by
// linking both sides into one binary.
//
// Unit conversion notes (see the project's own analysis before writing
// this): legacy's own XInteractor uses raw 1/r (see e.g. its R = 1 /
// e12.norm()), with no Coulomb-constant prefactor anywhere -- confirming
// legacy uses the same atomic-unit-style convention (e=1, 4*pi*eps0=1)
// as the new PolarSite/EwaldRealSpaceInteractor code, differing only in
// length scale (nm vs bohr). Position and induced dipole (charge*length)
// both convert the same way: value_bohr = value_nm * nm2bohr. Field
// (charge/length^2) converts the OPPOSITE direction, since it scales as
// the inverse square of length: value_bohr = value_nm * bohr2nm^2 (NOT
// nm2bohr^2 -- easy to get backwards, worth being explicit about).
//
// Sign convention note, confirmed by tracing XInteractor::BiasStat's 3-
// argument overload (the one actually called by PolarBackground) and
// XInteractor::FieldPerm_At_By together: BiasStat sets e12 = s +
// pol2.getPos() - pol1.getPos(), i.e. e12 points from site 1 toward site
// 2; FieldPerm_At_By then adds pol2's own multipole contribution into
// pol1.FPx/FPy/FPz proportional to +e12 (via T1x_00()=R2*rax with
// rax=e12(0)). The true physical field at site 1 from a positive charge
// at site 2 must point AWAY from site 2, i.e. in the -e12 direction --
// so APolarSite's own FPx/FPy/FPz stores the NEGATIVE of the physical
// field, the same "V=-E" convention independently confirmed elsewhere in
// this project for eeInteractor's own V() (see
// EwaldPeriodicDipoleOperator's own class documentation for that
// derivation). EwaldRealSpaceInteractor/EwaldRealSpaceSum/
// EwaldReciprocalSpaceSum's own V(), by contrast, is the genuine
// un-negated physical field (confirmed via test_ewaldrealspaceinteractor
// .cc's own explicit CalcStaticEnergy formula). The permanent field
// output here is therefore negated relative to getFieldP()'s own raw
// value, to match hdf5_dump's own (positive-field) convention -- this
// was caught the hard way, as a near-exactly-180-degree angle in
// compare_ewald.py's own output across essentially every site, only
// after two other hypotheses (quadrupoles, anisotropic polarizability)
// had already been ruled out by checking the actual force field. Do not
// revert this negation without rederiving the sign from the actual
// XInteractor code, not from intuition -- this exact mistake was made
// (and caught) more than once already in this project.
//
// Usage: ptop_dump <input.ptop> <output.csv>

// Standard includes
#include <fstream>
#include <iostream>

// Local VOTCA includes
#include "votca/tools/constants.h"
#include "votca/xtp/ewald/polartop.h"

using namespace votca;
using namespace votca::xtp;

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input.ptop> <output.csv>"
              << std::endl;
    return 1;
  }

  const double nm2bohr = tools::conv::nm2bohr;
  const double bohr2nm = tools::conv::bohr2nm;
  const double field_nm2bohr = bohr2nm * bohr2nm;  // see unit note above

  PolarTop ptop;
  ptop.LoadFromDrive(argv[1]);

  std::ofstream out(argv[2]);
  if (!out) {
    std::cerr << "Could not open output file: " << argv[2] << std::endl;
    return 1;
  }
  out << "segment_id,site_id,x_bohr,y_bohr,z_bohr,"
      << "FPx_au,FPy_au,FPz_au,U1x_au,U1y_au,U1z_au\n";

  std::size_t n_sites = 0;
  for (PolarSeg* seg : ptop.BGN()) {
    for (APolarSite* site : *seg) {
      vec pos = site->getPos() * nm2bohr;
      // Negated -- see the sign convention note above.
      vec fp = -site->getFieldP() * field_nm2bohr;
      vec u1 = site->getU1() * nm2bohr;
      out << seg->getId() << "," << site->getId() << "," << pos.x() << ","
          << pos.y() << "," << pos.z() << "," << fp.x() << "," << fp.y()
          << "," << fp.z() << "," << u1.x() << "," << u1.y() << ","
          << u1.z() << "\n";
      ++n_sites;
    }
  }

  std::cout << "Wrote " << n_sites << " sites from " << ptop.BGN().size()
            << " segments to " << argv[2] << std::endl;

  return 0;
}
