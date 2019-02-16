/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef __QMPACKAGEFACTORY__H
#define __QMPACKAGEFACTORY__H

#include <votca/tools/objectfactory.h>
#include <votca/xtp/qmpackage.h>

namespace votca {
namespace xtp {

class QMPackageFactory : public tools::ObjectFactory<std::string, QMPackage> {
 private:
  QMPackageFactory() {}

 public:
  static void RegisterAll(void);
  friend QMPackageFactory &QMPackages();
};

inline QMPackageFactory &QMPackages() {
  static QMPackageFactory _instance;
  return _instance;
}

}  // namespace xtp
}  // namespace votca

#endif /* __QMPACKAGEFACTORY__H */
