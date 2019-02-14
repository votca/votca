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

#ifndef _VOTCA_XTP_QMSANDBOX_H
#define _VOTCA_XTP_QMSANDBOX_H

#include <stdio.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmtool.h>

#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/qmpackagefactory.h>

namespace votca {
namespace xtp {

class QMSandbox : public QMTool {
 public:
  QMSandbox(){};
  ~QMSandbox(){};

  std::string Identify() { return "qmsandbox"; }

  void Initialize(tools::Property* options);
  bool Evaluate();

 private:
  std::string _orbfile;
  std::string _espfile;
  std::string _mpsfiled;
  std::string _mpsfileds;
  std::string _mpsfileq;
  std::string _mpsfileqs;
};

void QMSandbox::Initialize(tools::Property* options) {

  // update options with the VOTCASHARE defaults
  // UpdateWithDefaults( options, "xtp" );

  std::string key = "options." + Identify();

  _orbfile = options->ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".orbfile");
  _espfile = options->ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".espfile");
  _mpsfiled = options->ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".dipole");

  _mpsfileds = options->ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".dipole_split");
  _mpsfileq = options->ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".quadrupole");
  _mpsfileqs = options->ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".quadrupole_split");
}

bool QMSandbox::Evaluate() { return true; }

}  // namespace xtp
}  // namespace votca

#endif
