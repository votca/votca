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

#include <boost/program_options.hpp>
#include <iostream>
#include <votca/csg/version.h>
#include <votca/tools/akimaspline.h>
#include <votca/tools/cubicspline.h>
#include <votca/tools/linspline.h>
#include <votca/tools/spline.h>
#include <votca/tools/table.h>
#include <votca/tools/tokenizer.h>

using namespace std;
namespace po = boost::program_options;
using namespace votca::csg;
using namespace votca::tools;

void help_text() {
  votca::csg::HelpTextHeader("csg_resample");
  cout << "Change grid and interval of any sort of table files.\n"
          "Mainly called internally by inverse script, can also be\n"
          "used to manually prepare input files for coarse-grained\n"
          "simulations.\n\n";
}

void check_option(po::options_description &desc, po::variables_map &vm,
                  const string &option) {
  if (!vm.count(option)) {
    cout << "csg_resample \n\n";
    cout << desc << endl << "parameter " << option << " is not specified\n";
    exit(1);
  }
}

int main(int argc, char **argv) {

  string in_file, out_file, grid, fitgrid, comment, type, boundaries;
  Spline *spline = nullptr;
  Table in, out, der;
  // program options
  po::options_description desc("Allowed options");

  desc.add_options()("help", "produce this help message")(
      "in", po::value<string>(&in_file), "table to read")(
      "out", po::value<string>(&out_file), "table to write")(
      "derivative", po::value<string>(), "table to write")(
      "grid", po::value<string>(&grid),
      "new grid spacing (min:step:max). If 'grid' is specified only, "
      "interpolation is performed.")(
      "type", po::value<string>(&type)->default_value("akima"),
      "[cubic|akima|linear]. If option is not specified, the default type "
      "'akima' is assumed.")("fitgrid", po::value<string>(&fitgrid),
                             "specify fit grid (min:step:max). If 'grid' and "
                             "'fitgrid' are specified, a fit is performed.")(
      "nocut",
      "Option for fitgrid: Normally, values out of fitgrid boundaries are cut "
      "off. If they shouldn't, choose --nocut.")(
      "comment", po::value<string>(&comment),
      "store a comment in the output table")(
      "boundaries", po::value<string>(&boundaries),
      "(natural|periodic|derivativezero) sets boundary conditions");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (po::error &err) {
    cout << "error parsing command line: " << err.what() << endl;
    return -1;
  }

  try {
    // does the user want help?
    if (vm.count("help")) {
      help_text();
      cout << desc << endl;
      return 0;
    }

    check_option(desc, vm, "in");
    check_option(desc, vm, "out");

    if (!(vm.count("grid") || vm.count("fitgrid"))) {
      cout << "Need grid for interpolation or fitgrid for fit.\n";
      return 1;
    }

    if ((!vm.count("grid")) && vm.count("fitgrid")) {
      cout << "Need a grid for fitting as well.\n";
      return 1;
    }

    double min, max, step;
    {
      Tokenizer tok(grid, ":");
      vector<string> toks;
      tok.ToVector(toks);
      if (toks.size() != 3) {
        cout << "wrong range format, use min:step:max\n";
        return 1;
      }
      min = stod(toks[0]);
      step = stod(toks[1]);
      max = stod(toks[2]);
    }

    in.Load(in_file);

    if (vm.count("type")) {
      if (type == "cubic") {
        spline = new CubicSpline();
      } else if (type == "akima") {
        spline = new AkimaSpline();
      } else if (type == "linear") {
        spline = new LinSpline();
      } else {
        throw std::runtime_error("unknown type");
      }
    }
    spline->setBC(Spline::splineNormal);

    if (vm.count("boundaries")) {
      if (boundaries == "periodic") {
        spline->setBC(Spline::splinePeriodic);
      }
      if (boundaries == "derivativezero") {
        spline->setBC(Spline::splineDerivativeZero);
      }
      // default: normal
    }

    // in case fit is specified
    if (vm.count("fitgrid")) {
      Tokenizer tok(fitgrid, ":");
      vector<string> toks;
      tok.ToVector(toks);
      if (toks.size() != 3) {
        cout << "wrong range format in fitgrid, use min:step:max\n";
        return 1;
      }
      double sp_min, sp_max, sp_step;
      sp_min = stod(toks[0]);
      sp_step = stod(toks[1]);
      sp_max = stod(toks[2]);
      cout << "doing " << type << " fit " << sp_min << ":" << sp_step << ":"
           << sp_max << endl;

      // cut off any values out of fitgrid boundaries (exception: do nothing in
      // case of --nocut)
      Eigen::VectorXd x_copy;
      Eigen::VectorXd y_copy;
      if (!vm.count("nocut")) {
        // determine vector size
        int minindex = -1, maxindex = -1;
        for (int i = 0; i < in.x().size(); i++) {
          if (in.x(i) < sp_min) {
            minindex = i;
          }
          if (in.x(i) < sp_max) {
            maxindex = i;
          }
        }
        // copy data values in [sp_min,sp_max] into new vectors
        minindex++;
        x_copy = Eigen::VectorXd::Zero(maxindex - minindex + 1);
        y_copy = Eigen::VectorXd::Zero(maxindex - minindex + 1);
        for (int i = minindex; i <= maxindex; i++) {
          x_copy(i - minindex) = in.x(i);
          y_copy(i - minindex) = in.y(i);
        }
      }

      // fitting
      spline->GenerateGrid(sp_min, sp_max, sp_step);
      try {
        if (vm.count("nocut")) {
          spline->Fit(in.x(), in.y());
        } else {
          spline->Fit(x_copy, y_copy);
        }
      } catch (const char *message) {
        if (strcmp("qrsolve_zero_column_in_matrix", message)) {
          throw std::runtime_error(
              "error in Linalg::linalg_qrsolve : Not enough data for fit, "
              "please adjust grid (zero row in fit matrix)");
        } else if (strcmp("constrained_qrsolve_zero_column_in_matrix",
                          message)) {
          throw std::runtime_error(
              "error in Linalg::linalg_constrained_qrsolve : Not enough data "
              "for fit, please adjust grid (zero row in fit matrix)");
        } else
          throw std::runtime_error(
              "Unknown error in csg_resample while fitting.");
      }
    } else {
      // otherwise do interpolation (default = cubic)
      try {
        spline->Interpolate(in.x(), in.y());
      } catch (const char *message) {
        if (strcmp("qrsolve_zero_column_in_matrix", message)) {
          throw std::runtime_error(
              "error in Linalg::linalg_qrsolve : Not enough data, please "
              "adjust grid (zero row in fit matrix)");
        } else if (strcmp("constrained_qrsolve_zero_column_in_matrix",
                          message)) {
          throw std::runtime_error(
              "error in Linalg::linalg_constrained_qrsolve : Not enough data, "
              "please adjust grid (zero row in fit matrix)");
        } else
          throw std::runtime_error(
              "Unknown error in csg_resample while interpolating.");
      }
    }

    out.GenerateGridSpacing(min, max, step);
    spline->Calculate(out.x(), out.y());

    // store a comment line
    if (vm.count("comment")) {
      out.set_comment(comment);
    }
    out.y() = out.y();
    out.flags() = std::vector<char>(out.flags().size(), 'o');

    der.GenerateGridSpacing(min, max, step);
    der.flags() = std::vector<char>(der.flags().size(), 'o');

    int i = 0;
    for (i = 0; out.x(i) < in.x(0) && i < out.size(); ++i)
      ;

    int j = 0;
    for (; i < out.size(); ++i) {
      for (; j < in.size(); ++j)
        if (in.x(j) >= out.x(i) ||
            fabs(in.x(j) - out.x(i)) < 1e-12)  // fix for precison errors
          break;
      if (in.size() == j) break;
      out.flags(i) = in.flags(j);
      der.flags(i) = in.flags(j);
    }

    out.Save(out_file);

    if (vm.count("derivative")) {
      spline->CalculateDerivative(der.x(), der.y());

      der.Save(vm["derivative"].as<string>());
    }

    delete spline;
  } catch (std::exception &error) {
    cerr << "an error occurred:\n" << error.what() << endl;
    return -1;
  }
  return 0;
}
