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

#ifndef __VOTCA_XTP_QMDATABASE2_H
#define __VOTCA_XTP_QMDATABASE2_H

#include <votca/tools/database.h>

using namespace votca::tools;

namespace votca {
namespace xtp {

/**
 * \brief the state database
 *
 * This class contains management of state databases. It creates
 * new databases and opens or upgrades existing ones.
 */
class QMDatabase : public Database {
 public:
  /**
   * \brief Create the database scheme
   *
   * This function is called when a database is created.
   * All tables and triggers should be added here.
   */
  void onCreate();
};

}  // namespace xtp
}  // namespace votca

#endif /* QMDATABASE2_H */
