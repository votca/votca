/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include <string>
#include <vector>
#include <votca/tools/eigen.h>
#include <votca/tools/types.h>

namespace votca {
namespace tools {

/**
    \brief class to store tables like rdfs, tabulated potentials, etc

    \think about weather to make this a template, can be used in histogram
    as well, of for counting with integers...

 */
class Table {
 public:
  void clear();

  void GenerateGridSpacing(double min, double max, double spacing);
  void resize(Index N);
  Index size() const { return _x.size(); }

  double &x(Index i) { return _x[i]; }
  double &y(Index i) { return _y[i]; }

  const double &x(Index i) const { return _x[i]; }
  const double &y(Index i) const { return _y[i]; }
  char &flags(Index i) { return _flags[i]; }
  double &yerr(Index i) { return _yerr[i]; }

  void set(const Index &i, const double &x, const double &y) {
    _x[i] = x;
    _y[i] = y;
  }
  void set(const Index &i, const double &x, const double &y,
           const char &flags) {
    _x[i] = x;
    _y[i] = y;
    _flags[i] = flags;
  }
  void set(const Index &i, const double &x, const double &y, const char &flags,
           const double &yerr) {
    _x[i] = x;
    _y[i] = y;
    _flags[i] = flags;
    _yerr[i] = yerr;
  }

  void set_comment(const std::string comment) {
    _has_comment = true;
    _comment_line = comment;
  }

  void Load(std::string filename);
  void Save(std::string filename) const;

  void Smooth(Index Nsmooth);

  bool GetHasYErr() { return _has_yerr; }
  void SetHasYErr(bool has_yerr) { _has_yerr = has_yerr; }

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

  Eigen::VectorXd &x() { return _x; }
  Eigen::VectorXd &y() { return _y; }
  std::vector<char> &flags() { return _flags; }
  Eigen::VectorXd &yerr() { return _yerr; }

  void push_back(double x, double y, char flags = ' ');

  const std::string &getErrorDetails() { return _error_details; }

  void setErrorDetails(std::string str) { _error_details = str; }

 private:
  Eigen::VectorXd _x;
  Eigen::VectorXd _y;
  std::vector<char> _flags;
  Eigen::VectorXd _yerr;
  std::string _error_details = "";

  bool _has_yerr = false;
  bool _has_comment = false;

  friend std::ostream &operator<<(std::ostream &out, const Table &v);
  friend std::istream &operator>>(std::istream &in, Table &t);

  std::string _comment_line;
};

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_TABLE_H
