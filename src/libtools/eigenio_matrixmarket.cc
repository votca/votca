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

// Local VOTCA includes
#include "votca/tools/eigen.h"
#include "votca/tools/tokenizer.h"
#include "votca/tools/types.h"

namespace votca {
namespace tools {

namespace EigenIO_MatrixMarket {

Eigen::VectorXd ReadVector(const std::string& filename) {

  Eigen::VectorXd output;

  bool success = Eigen::loadMarketVector(output, filename);
  if (!success) {
    throw std::runtime_error("Loading Vector from " + filename + " failed");
  }
  return output;
}

void WriteVector(const std::string& filename, const Eigen::VectorXd& output) {
  bool success = Eigen::saveMarketVector(output, filename);
  if (!success) {
    throw std::runtime_error("Writing Vector to " + filename + " failed");
  }
}

void WriteMatrix(const std::string& filename, const Eigen::MatrixXd& output) {

  std::ofstream ofs;
  ofs.open(filename, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Could not create " + filename);
  }

  ofs << "%%MatrixMarket matrix array real general\n";
  ofs << output.rows() << " " << output.cols() << "\n";
  const Eigen::Map<const Eigen::VectorXd> temp(output.data(), output.size());
  ofs << temp;
  ofs.close();
}

Eigen::MatrixXd ReadMatrix(const std::string& filename) {
  int sym;
  bool iscomplex;
  bool isvector;
  bool success = Eigen::getMarketHeader(filename, sym, iscomplex, isvector);
  if (!success) {
    throw std::runtime_error("Could not read " + filename);
  }
  if (sym != 0) {
    throw std::runtime_error("Only supports reading in general matrices");
  }
  if (iscomplex) {
    throw std::runtime_error("Only supports reading in real matrices");
  }

  if (!isvector) {
    throw std::runtime_error(
        "Use the eigen method `loadMarket` for  sparse data");
  }

  std::ifstream in(filename, std::ios::in);
  if (!in) {
    throw std::runtime_error("Could not open " + filename);
  }

  std::string line;
  Index cols = 0;
  Index rows = 0;
  do {  // Skip comments
    std::getline(in, line);
    assert(in.good());
  } while (line[0] == '%');
  std::istringstream newline(line);
  newline >> rows >> cols;
  assert(rows > 0 && cols > 0);

  Index i = 0;
  Index n = rows * cols;
  std::vector<double> entries;
  entries.reserve(n);
  double value;
  while (std::getline(in, line) && (i < n)) {
    std::istringstream newline2(line);
    newline2 >> value;
    entries.push_back(value);
    i++;
  }
  in.close();

  return Eigen::Map<Eigen::MatrixXd>(entries.data(), rows, cols);
}

}  // namespace EigenIO_MatrixMarket

}  // namespace tools
}  // namespace votca
