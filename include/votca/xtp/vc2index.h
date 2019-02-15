/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_VC2INDEX_H
#define _VOTCA_XTP_VC2INDEX_H

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
  vc2index(int vmin, int cmin, int ctotal)
      : _vmin(vmin), _cmin(cmin), _ctotal(ctotal){};

  inline int I(int v, int c) const { return _ctotal * (v - _vmin) + c - _cmin; }

  inline int v(int index) const { return index / _ctotal + _vmin; }

  inline int c(int index) const { return index % _ctotal + _cmin; }

 private:
  int _vmin;
  int _cmin;
  int _ctotal;
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_VC2INDEX_H */
