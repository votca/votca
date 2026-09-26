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
// manpage, or sphinx-doc infrastructure). Dumps an EwaldRegistry HDF5
// checkpoint (as written by EwaldBackground) as a CSV, one row per site,
// one row per site. The schema was originally shared with ptop_dump, the
// matching dumper for the legacy xtp/ewald/ PolarTop, so that the two
// could be diffed directly; that tool went with the code it read, and
// this one remains as the way to inspect a background checkpoint.
//
// Positions, permanent field, and induced dipoles are already in this
// tool's own native units (bohr, atomic-unit-style) -- no conversion
// needed, unlike the nm->bohr/field conversion the legacy dumper did.
//
// The permanent-field columns (FPx/FPy/FPz) read out PolarSite's own V()
// accumulator. This only holds anything meaningful because
// EwaldBackground::Evaluate() now explicitly writes the permanent field
// (b) back into every site's V() right before the checkpoint write --
// neither of its own two code paths (the PCG solve, or the
// induce=false skip) otherwise leaves V() holding it at that point (see
// EwaldBackground's own comment at that write for the full trace). If
// you're reading this from a checkpoint written by an older version of
// EwaldBackground that predates that write-back, these columns will
// silently be zero instead of a real error -- worth checking the
// EwaldBackground version/log if these numbers look suspiciously flat.
//
// Usage: hdf5_dump <input.hdf5> <output.csv>

// Standard includes
#include <fstream>
#include <iostream>

// Local VOTCA includes
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/ewaldregistry.h"

using namespace votca;
using namespace votca::xtp;

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input.hdf5> <output.csv>"
              << std::endl;
    return 1;
  }

  CheckpointFile cpf(argv[1], CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader();

  EwaldRegistry registry;
  registry.ReadFromCpt(r);

  std::ofstream out(argv[2]);
  if (!out) {
    std::cerr << "Could not open output file: " << argv[2] << std::endl;
    return 1;
  }
  out << "segment_id,site_id,x_bohr,y_bohr,z_bohr,"
      << "FPx_au,FPy_au,FPz_au,U1x_au,U1y_au,U1z_au\n";

  std::size_t n_sites = 0;
  std::size_t n_segments = 0;
  for (Index id : registry.AllIds()) {
    if (!registry.Has(id, EwaldChargeState::Neutral)) {
      continue;
    }
    const PolarSegment& segment = registry.Get(id, EwaldChargeState::Neutral);
    ++n_segments;
    for (const PolarSite& site : segment) {
      const Eigen::Vector3d& pos = site.getPos();
      const Eigen::Vector3d& fp = site.V();
      const Eigen::Vector3d& u1 = site.getInducedDipole();
      out << id << "," << site.getId() << "," << pos.x() << "," << pos.y()
          << "," << pos.z() << "," << fp.x() << "," << fp.y() << "," << fp.z()
          << "," << u1.x() << "," << u1.y() << "," << u1.z() << "\n";
      ++n_sites;
    }
  }

  std::cout << "Wrote " << n_sites << " sites from " << n_segments
            << " segments to " << argv[2] << std::endl;

  return 0;
}
