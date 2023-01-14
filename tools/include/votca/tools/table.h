/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef VOTCA_TOOLS_TABLE_H
#define VOTCA_TOOLS_TABLE_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "types.h"

namespace votca {
namespace tools {

/**
    \brief class to store tables like rdfs, tabulated potentials, etc

 */
class Table {
 public:
  void clear();

  void GenerateGridSpacing(double min, double max, double spacing);
  void resize(Index N);
  Index size() const { return x_.size(); }

  double &x(Index i) { return x_[i]; }
  double &y(Index i) { return y_[i]; }

  const double &x(Index i) const { return x_[i]; }
  const double &y(Index i) const { return y_[i]; }
  char &flags(Index i) { return flags_[i]; }
  double &yerr(Index i) { return yerr_[i]; }

  void set(const Index &i, const double &x, const double &y) {
    x_[i] = x;
    y_[i] = y;
  }
  void set(const Index &i, const double &x, const double &y,
           const char &flags) {
    x_[i] = x;
    y_[i] = y;
    flags_[i] = flags;
  }
  void set(const Index &i, const double &x, const double &y, const char &flags,
           const double &yerr) {
    x_[i] = x;
    y_[i] = y;
    flags_[i] = flags;
    yerr_[i] = yerr;
  }

  void set_comment(const std::string comment) {
    has_comment_ = true;
    comment_line_ = comment;
  }

  void Load(std::string filename);
  void Save(std::string filename) const;

  void Smooth(Index Nsmooth);

  bool GetHasYErr() { return has_yerr_; }
  void SetHasYErr(bool has_yerr) { has_yerr_ = has_yerr; }

  /**
   * \brief Gets the maximum value in the y column
   * \return - max value
   */
  double getMaxY() const;
  /**
   * \brief Gets the minimum value in the y column
   * \return - min value
   */
  double getMinY() const;
  /**
   * \brief Gets the maximum value in the x column
   * \return - max value
   */
  double getMaxX() const;
  /**
   * \brief Gets the minimum value in the x column
   * \return - min value
   */
  double getMinX() const;

  Eigen::VectorXd &x() { return x_; }
  Eigen::VectorXd &y() { return y_; }
  std::vector<char> &flags() { return flags_; }
  Eigen::VectorXd &yerr() { return yerr_; }

  void push_back(double x, double y, char flags = ' ');

  const std::string &getErrorDetails() { return error_details_; }

  void setErrorDetails(std::string str) { error_details_ = str; }

 private:
  Eigen::VectorXd x_;
  Eigen::VectorXd y_;
  std::vector<char> flags_;
  Eigen::VectorXd yerr_;
  std::string error_details_ = "";

  bool has_yerr_ = false;
  bool has_comment_ = false;

  friend std::ostream &operator<<(std::ostream &out, const Table &t);
  friend std::istream &operator>>(std::istream &in, Table &t);

  std::string comment_line_;
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_TABLE_H
