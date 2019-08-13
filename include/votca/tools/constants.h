/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_TOOLS_CONSTANTS_H
#define VOTCA_TOOLS_CONSTANTS_H

#include <boost/math/constants/constants.hpp>
#include <cmath>

namespace votca {
namespace tools {

namespace conv {

// mathematical constants
const double Pi = boost::math::constants::pi<double>();

// natural constants
// Boltzmann Factor eV/K
const double kB = 8.617332478E-5;
// Planck's Constant eV*s
const double hbar = 6.5821192815E-16;

// length conversions
// votca xtp-uses for any conversions the following scheme unitA2unitB
const double bohr2nm = 0.052917721092;
const double nm2bohr = 18.897259886;
const double ang2bohr = 1.8897259886;
const double bohr2ang = 1.0 / 1.8897259886;
const double nm2ang = 10.0;
const double ang2nm = 0.1;

const double hrt2ev = 27.21138602;
const double ev2hrt = 1.0 / 27.21138602;

// 1 eV = 96.485 Kj/mol
const double ev2kj_per_mol = 96.485;

}  // namespace conv

namespace topology_constants {

/// Used to indicate that a valid element variable has not been assigned
const std::string unassigned_element = "unassigned";

/// Used to indicate that a valid bead type variable has not been assigned
const std::string unassigned_bead_type = "unassigned";

/// Used to indicate that a valid residue type variable has not been assigned
const std::string unassigned_residue_type = "unassigned";

/// Used to indicate that a valid molecule type variable has not been assigned
const std::string unassigned_molecule_type = "unassigned";

/// Used to indicate that a valid segment type variable has not been assigned
const std::string unassigned_segment_type = "unassigned";

const std::string unassigned_atom_container_type = "unassigned";

const int unassigned_atom_container_id = -1;
/// Used to indicate a valid residue id has not been assigned
const int unassigned_residue_id = -1;

/// Used to indicate a valid molecule id has not been assigned
const int unassigned_molecule_id = -1;

/// Used to indicate a valid molecule id has not been assigned
const int unassigned_bead_id = -1;

/// Used to indicate a valid segment id has not been assigned
const int unassigned_segment_id = -1;
}  // namespace topology_constants
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_CONSTANTS_H
