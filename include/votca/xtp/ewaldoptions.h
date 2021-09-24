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

#pragma once
#ifndef VOTCA_XTP_EWALDOPTIONS_H
#define VOTCA_XTP_EWALDOPTIONS_H
enum class Shape { xyslab, cube, sphere };

struct EwaldOptions {
  double k_cutoff;
  double r_cutoff;
  double alpha;
  double sharpness;
  Shape shape;

  bool operator==(const EwaldOptions& other){
    bool equal = true;
    equal = equal && (this->alpha == other.alpha);
    equal = equal && (this->k_cutoff == other.k_cutoff);
    equal = equal && (this->r_cutoff == other.r_cutoff);
    equal = equal && (this->sharpness == other.sharpness);
    equal = equal && (this->shape == other.shape);
    return equal;
  }

  bool operator!=(const EwaldOptions& other){
    return !operator==(other);
  }

};

#endif