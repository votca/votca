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

/*
 *  csg_imcrepack repacks the matri for the imc update. It either removes zero
 * line entries (pack) or extracts the dU from the big file (unpack)
 *
 */

#include <boost/program_options.hpp>
#include <iostream>
#include <votca/csg/imcio.h>
#include <votca/csg/version.h>
#include <votca/tools/table.h>
#include <votca/tools/tokenizer.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

void help_text(void) {
  votca::csg::HelpTextHeader("csg_imcrepack");
  cout << "This program is internally called by inversion scripts to kick out\n"
          "zero entries in matrix for inverse Monte Carlo. It also extracts "
          "the\n"
          "single potential updates out of the full solution.\n\n";
}

using namespace std;
namespace po = boost::program_options;

void check_option(po::options_description &desc, po::variables_map &vm,
                  const string &option) {
  if (!vm.count(option)) {
    cout << "csg_imcrepack \n\n";
    cout << desc << endl << "parameter " << option << " is not specified\n";
    exit(1);
  }
}

int main(int argc, char **argv) {
  string name_in, name_out, name_unpack;
  // program options
  po::options_description desc("Allowed options");

  desc.add_options()("in", po::value<string>(&name_in), "files to read")(
      "out", po::value<string>(&name_out), "files to write")(
      "unpack", po::value<string>(&name_unpack),
      "extract all tables from this file")("help", "display help message");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (po::error &err) {
    cout << "error parsing command line: " << err.what() << endl;
    return -1;
  }

  // does the user want help?
  if (vm.count("help")) {
    help_text();
    cout << desc << endl;
    return 0;
  }

  check_option(desc, vm, "in");

  std::vector<std::pair<std::string, RangeParser> > ranges =
      imcio_read_index(name_in + ".idx");

  if (vm.count("unpack")) {
    Table tbl_in;
    tbl_in.Load(name_unpack);
    for (std::pair<std::string, RangeParser> &range : ranges) {
      Table tbl;
      for (int r : range.second) {
        tbl.push_back(tbl_in.x(r - 1), tbl_in.y(r - 1), 'i');
      }
      tbl.Save(range.first + ".dpot.imc");
    }
  } else {
    check_option(desc, vm, "out");

    int beg = 1;
    int end = 0;
    std::list<int> list;
    Table dS;
    dS.Load(name_in + ".imc");
    Eigen::MatrixXd gmc = imcio_read_matrix(name_in + ".gmc");

    for (std::pair<std::string, RangeParser> &range : ranges) {

      for (int r : range.second) {
        for (int i = 0; i < gmc.rows(); ++i) {
          if (((gmc.row(i).cwiseAbs().array()) > 1e-8).any()) {
            list.push_back(r - 1);
            end++;
          }
        }
      }
      RangeParser new_rp;
      new_rp.Add(beg, end);
      beg = end + 1;
      range.second = new_rp;
    }
    imcio_write_dS(name_out + ".imc", dS, &list);
    imcio_write_matrix(name_out + ".gmc", gmc, &list);
    imcio_write_index(name_out + ".idx", ranges);
  }
}
