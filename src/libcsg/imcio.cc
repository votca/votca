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

#include "../../include/votca/csg/imcio.h"
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <votca/tools/getline.h>
#include <votca/tools/rangeparser.h>
#include <votca/tools/table.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

using group_matrix = Eigen::MatrixXd;
using namespace std;

void imcio_write_dS(const string &file, const tools::Table &dS,
                    const std::list<Index> *list) {
  // write the dS
  ofstream out_dS;
  out_dS.open(file);
  out_dS << setprecision(8);
  if (!out_dS) {
    throw runtime_error(string("error, cannot open file ") + file);
  }

  if (list == nullptr) {
    for (Index i = 0; i < dS.size(); ++i) {
      out_dS << dS.x(i) << " " << dS.y(i) << endl;
    }
  } else {
    for (Index i : *list) {
      out_dS << dS.x(i) << " " << dS.y(i) << endl;
    }
  }

  out_dS.close();
  cout << "written " << file << endl;
}

void imcio_write_matrix(const string &file, const Eigen::MatrixXd &gmc,
                        const std::list<Index> *list) {
  ofstream out_A;
  out_A.open(file);
  out_A << setprecision(8);

  if (!out_A) {
    throw runtime_error(string("error, cannot open file ") + file);
  }

  if (list == nullptr) {
    for (Index i = 0; i < gmc.rows(); ++i) {
      for (Index j = 0; j < gmc.cols(); ++j) {
        out_A << gmc(i, j) << " ";
      }
      out_A << endl;
    }
  } else {
    for (Index i : *list) {
      for (Index j : *list) {
        out_A << gmc(i, j) << " ";
      }
      out_A << endl;
    }
  }
  out_A.close();
  cout << "written " << file << endl;
}

void imcio_write_index(
    const string &file,
    const std::vector<std::pair<std::string, tools::RangeParser> > &ranges) {
  // write the index
  ofstream out_idx;
  out_idx.open(file);

  if (!out_idx) {
    throw runtime_error(string("error, cannot open file ") + file);
  }

  for (const auto &range : ranges) {
    out_idx << range.first << " " << range.second << endl;
  }
  out_idx.close();
  cout << "written " << file << endl;
}

Eigen::MatrixXd imcio_read_matrix(const std::string &filename) {
  std::ifstream intt;
  intt.open(filename);
  if (!intt) {
    throw std::runtime_error(std::string("error, cannot open file ") +
                             filename);
  }

  std::string line;
  std::vector<double> result;
  Index numrows = 0;
  size_t numcols = 0;
  while (getline(intt, line)) {
    if (line[0] == '#') {
      continue;
    }
    votca::tools::Tokenizer tok(line, " \t");
    std::vector<std::string> tokens = tok.ToVector();
    if (numrows == 0) {
      numcols = tokens.size();
    } else if (numcols != tokens.size()) {
      throw std::runtime_error(
          "Matrix has not the same number of entries in each row.");
    }
    for (const std::string &s : tokens) {
      result.push_back(std::stod(s));
    }
    numrows++;
  }
  intt.close();

  return Eigen::Map<Eigen::MatrixXd>(result.data(), numrows, numcols);
}

std::vector<std::pair<std::string, tools::RangeParser> > imcio_read_index(
    const string &filename) {
  ifstream in;
  in.open(filename);
  if (!in) {
    throw runtime_error(string("error, cannot open file ") + filename);
  }

  std::vector<std::pair<std::string, tools::RangeParser> > indeces;
  string line;
  // read till the first data line
  while (getline(in, line)) {
    // remove comments and xmgrace stuff
    line = line.substr(0, line.find("#"));
    line = line.substr(0, line.find("@"));

    boost::trim(line);
    size_t found = line.find(" ");
    if (found == string::npos) {
      throw runtime_error(string("wrong format in ") + filename);
    }

    string name = line.substr(0, found);
    string range = line.substr(found);
    tools::RangeParser rp;
    rp.Parse(range);
    indeces.push_back(std::pair<std::string, tools::RangeParser>(name, rp));
  }
  in.close();
  return indeces;
}

}  // namespace csg
}  // namespace votca
