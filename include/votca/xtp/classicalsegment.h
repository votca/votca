/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef VOTCA_XTP_CLASSICALSEGMENT_H
#define VOTCA_XTP_CLASSICALSEGMENT_H

#include <votca/xtp/atomcontainer.h>
#include <votca/xtp/polarsite.h>
#include <votca/xtp/staticsite.h>

namespace votca {
namespace xtp {
template <class T>
class ClassicalSegment : public AtomContainer<T> {
 public:
  ClassicalSegment(std::string name, int id) : AtomContainer<T>(name, id){};

  ClassicalSegment(CheckpointReader& r) : AtomContainer<T>(r){};

  void LoadFromFile(std::string filename);

  void WriteMPS(std::string filename, std::string header) const;

  double CalcTotalQ() const;

  std::string identify() const;

  Eigen::Vector3d CalcDipole() const;

  friend std::ostream& operator<<(std::ostream& out,
                                  const ClassicalSegment<T>& container) {
    out << container.getId() << " " << container.getName() << "\n";
    for (const T& atom : container) {
      out << atom;
    }
    out << std::endl;
    return out;
  }
};

typedef ClassicalSegment<PolarSite> PolarSegment;
typedef ClassicalSegment<StaticSite> StaticSegment;

}  // namespace xtp
}  // namespace votca

#endif /* VOTCA_XTP_CLASSICALSEGMENT_H */
