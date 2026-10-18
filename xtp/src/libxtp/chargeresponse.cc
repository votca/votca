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
#include <fstream>
#include <sstream>
#include <stdexcept>

// Third party includes
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/getline.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/chargeresponse.h"

namespace votca {
namespace xtp {

namespace {
// Relative tolerance of the checks of a .crm file: the fitting script
// writes 9 significant digits, so a correct file is symmetric and
// charge conserving to about 1e-8 of its largest element.
constexpr double crm_tolerance = 1e-6;
}  // namespace

ChargeResponse::ChargeResponse(std::vector<std::string> elements,
                               Eigen::MatrixXd pi)
    : elements_(std::move(elements)), pi_(std::move(pi)) {
  Check("ChargeResponse constructor");
}

void ChargeResponse::LoadFromFile(const std::string& filename) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    throw std::runtime_error("File:" + filename + " could not be opened");
  }
  // all non-comment tokens in file order
  std::vector<std::string> tokens;
  std::string line;
  while (in.good()) {
    tools::getline(in, line);
    tools::Tokenizer toker(line, " \t");
    std::vector<std::string> split = toker.ToVector();
    if (split.empty() || split[0].substr(0, 1) == "!") {
      continue;
    }
    tokens.insert(tokens.end(), split.begin(), split.end());
  }
  auto fail = [&](const std::string& what) {
    throw std::runtime_error("Charge-response file " + filename + ": " + what);
  };
  if (tokens.empty()) {
    fail("no data");
  }
  Index n = 0;
  try {
    n = boost::lexical_cast<Index>(tokens[0]);
  } catch (const boost::bad_lexical_cast&) {
    fail("the first entry must be the number of atoms, got '" + tokens[0] +
         "'");
  }
  if (n <= 0) {
    fail("number of atoms must be positive");
  }
  const std::size_t expected = 1 + 4 * std::size_t(n) + std::size_t(n * n);
  if (tokens.size() != expected) {
    fail("expected " + std::to_string(expected) + " entries for " +
         std::to_string(n) +
         " atoms (N, N lines 'El x y z', N x N matrix), "
         "found " +
         std::to_string(tokens.size()));
  }
  elements_.clear();
  std::size_t k = 1;
  try {
    for (Index i = 0; i < n; ++i) {
      elements_.push_back(tokens[k]);
      // the geometry of the fit is informational: Pi is rotation invariant
      // and the sites come from the .mps file
      for (Index c = 0; c < 3; ++c) {
        (void)boost::lexical_cast<double>(tokens[k + 1 + std::size_t(c)]);
      }
      k += 4;
    }
    pi_.resize(n, n);
    for (Index i = 0; i < n; ++i) {
      for (Index j = 0; j < n; ++j) {
        pi_(i, j) = boost::lexical_cast<double>(tokens[k++]);
      }
    }
  } catch (const boost::bad_lexical_cast&) {
    fail("could not read entry " + std::to_string(k) + " ('" + tokens[k] +
         "') as a number");
  }
  Check(filename);
}

void ChargeResponse::Check(const std::string& source) {
  auto fail = [&](const std::string& what) {
    throw std::runtime_error("Charge-response matrix from " + source + ": " +
                             what);
  };
  const Index n = Index(pi_.rows());
  if (pi_.cols() != n || n == 0) {
    fail("must be a non-empty square matrix");
  }
  if (Index(elements_.size()) != n) {
    fail("has " + std::to_string(n) + " rows but " +
         std::to_string(elements_.size()) + " atoms");
  }
  const double scale = pi_.cwiseAbs().maxCoeff();
  if (!(scale > 0.0) || !pi_.allFinite()) {
    fail("must be finite and not zero");
  }
  const double asym = (pi_ - pi_.transpose()).cwiseAbs().maxCoeff();
  if (asym > crm_tolerance * scale) {
    fail("is not symmetric (largest |Pi_ij - Pi_ji| " + std::to_string(asym) +
         ")");
  }
  const double rowsum = pi_.rowwise().sum().cwiseAbs().maxCoeff();
  if (rowsum > crm_tolerance * scale * double(n)) {
    fail("does not conserve charge (largest row sum " + std::to_string(rowsum) +
         ")");
  }
  // remove the rounding of the file: symmetrize, project out the constant
  // vector (charge conservation)
  const Eigen::MatrixXd P = Eigen::MatrixXd::Identity(n, n) -
                            Eigen::MatrixXd::Constant(n, n, 1.0 / double(n));
  pi_ = P * (0.5 * (pi_ + pi_.transpose())) * P;
  pi_ = 0.5 * (pi_ + pi_.transpose());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(pi_,
                                                    Eigen::EigenvaluesOnly);
  const double lowest = es.eigenvalues()(0);
  const double highest = es.eigenvalues()(n - 1);
  if (lowest < -crm_tolerance * highest) {
    fail("is not positive semidefinite (lowest eigenvalue " +
         std::to_string(lowest) + ", highest " + std::to_string(highest) + ")");
  }
}

ChargeResponse ChargeResponse::Permuted(const std::vector<Index>& order) const {
  const Index n = size();
  if (Index(order.size()) != n) {
    throw std::runtime_error(
        "ChargeResponse::Permuted: the charge response of " +
        std::to_string(n) + " atoms needs all of them, got " +
        std::to_string(order.size()) +
        " (a subset of a molecule does not conserve charge)");
  }
  std::vector<bool> seen(std::size_t(n), false);
  for (Index i : order) {
    if (i < 0 || i >= n || seen[std::size_t(i)]) {
      throw std::runtime_error(
          "ChargeResponse::Permuted: the atom order is not a permutation of "
          "0.." +
          std::to_string(n - 1));
    }
    seen[std::size_t(i)] = true;
  }
  ChargeResponse result;
  result.pi_.resize(n, n);
  for (Index i = 0; i < n; ++i) {
    result.elements_.push_back(elements_[std::size_t(order[std::size_t(i)])]);
    for (Index j = 0; j < n; ++j) {
      result.pi_(i, j) = pi_(order[std::size_t(i)], order[std::size_t(j)]);
    }
  }
  return result;
}

void ChargeResponse::WriteToCpt(CheckpointWriter& w) const {
  std::string elements;
  for (const std::string& e : elements_) {
    elements += (elements.empty() ? "" : " ") + e;
  }
  w(elements, "elements");
  w(pi_, "Pi");
}

void ChargeResponse::ReadFromCpt(CheckpointReader& r) {
  std::string elements;
  r(elements, "elements");
  std::istringstream in(elements);
  elements_.clear();
  std::string e;
  while (in >> e) {
    elements_.push_back(e);
  }
  r(pi_, "Pi");
  Check("checkpoint");
}

}  // namespace xtp
}  // namespace votca
