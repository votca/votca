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
#include <stdexcept>
#include <string>

// Local VOTCA includes
#include "votca/xtp/ewaldparameters.h"

namespace votca {
namespace xtp {

namespace {
// The shape travels as a string rather than as the enum's integer value.
// An integer would silently change meaning if an enumerator were ever
// inserted ahead of another, and a checkpoint written by an older build
// would then be read back as a different shape with no error -- the
// kind of mismatch that produces a plausible wrong answer.
std::string ShapeToString(EwaldShape shape) {
  switch (shape) {
    case EwaldShape::Cube:
      return "cube";
    case EwaldShape::Slab:
      return "slab";
  }
  throw std::runtime_error("EwaldParameters: unhandled EwaldShape");
}

EwaldShape ShapeFromString(const std::string& name) {
  if (name == "cube") {
    return EwaldShape::Cube;
  }
  if (name == "slab") {
    return EwaldShape::Slab;
  }
  throw std::runtime_error("EwaldParameters: unknown shape '" + name +
                           "' in checkpoint");
}
}  // namespace

void EwaldParameters::WriteToCpt(CheckpointWriter& w) const {
  w(alpha, "alpha");
  w(k_max, "k_max");
  w(r_min, "r_min");
  w(field_tol, "field_tol");
  w(thole_a, "thole_a");
  w(screening_factor, "screening_factor");
  w(ShapeToString(shape), "shape");
  w(box, "box");
}

void EwaldParameters::ReadFromCpt(CheckpointReader& r) {
  r(alpha, "alpha");
  r(k_max, "k_max");
  r(r_min, "r_min");
  r(field_tol, "field_tol");
  r(thole_a, "thole_a");
  r(screening_factor, "screening_factor");
  std::string shape_name;
  r(shape_name, "shape");
  shape = ShapeFromString(shape_name);
  r(box, "box");

  if (!(alpha > 0.0) || !(k_max > 0.0) || !(box.determinant() > 0.0)) {
    throw std::runtime_error(
        "EwaldParameters: checkpoint holds implausible Ewald parameters "
        "(alpha, k_max and the box volume must all be positive). The "
        "checkpoint is either from an incompatible version or corrupt.");
  }
}

}  // namespace xtp
}  // namespace votca
