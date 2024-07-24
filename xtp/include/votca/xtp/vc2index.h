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
#ifndef VOTCA_XTP_VC2INDEX_H
#define VOTCA_XTP_VC2INDEX_H

// VOTCA includes
#include <votca/tools/types.h>

namespace votca {
namespace xtp {
/**
 * \brief This class transforms a pair of indices v,c into a compound index I,
 * via I=ctotal*v+c the fast dimension is c. If you have a choice iterate over c
 * and v not over I.
 *
 *
 */
class vc2index {

 public:
  vc2index(Index vmin, Index cmin, Index ctotal)
      : vmin_(vmin), cmin_(cmin), ctotal_(ctotal) {};

  inline Index I(Index v, Index c) const {
    return ctotal_ * (v - vmin_) + (c - cmin_);
  }

  inline Index v(Index index) const { return (index / ctotal_ + vmin_); }

  inline Index c(Index index) const { return (index % ctotal_ + cmin_); }

 private:
  Index vmin_;
  Index cmin_;
  Index ctotal_;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_VC2INDEX_H
