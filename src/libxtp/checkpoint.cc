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

#include "gitversion.h"
#include <iostream>
#include <stdexcept>
#include <votca/xtp/checkpoint.h>

namespace votca {
namespace xtp {

using namespace hdf5_utils;

CheckpointFile::CheckpointFile(std::string fN)
    : _fileName(fN), _version(gitversion) {

  try {
    H5::Exception::dontPrint();
    _fileHandle = H5::H5File(_fileName, H5F_ACC_TRUNC);

    WriteAttribute(_fileHandle, gitversion, "Version");

    _child1 = _fileHandle.createGroup("/Child1");

    Eigen::Matrix<double, 10, 10> b;

    // WriteData(_child1, b, "b");

    Eigen::VectorXd c(10);

    WriteData(_child1, c, "c");

  } catch (H5::Exception error) {
    error.printError();
    throw std::runtime_error(error.getDetailMsg());
  }
};

std::string CheckpointFile::getFileName() { return _fileName; };
std::string CheckpointFile::getVersion() { return _version; };

}  // namespace xtp
}  // namespace votca
