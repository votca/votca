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

// Standard includes
#include <cstdio>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/getline.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/cubefile_writer.h"
#include "votca/xtp/orbitals.h"

// Local private VOTCA includes
#include "gencube.h"

namespace votca {
namespace xtp {

void GenCube::ParseOptions(const tools::Property& options) {

  orbfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".input", job_name_ + ".orb");
  output_file_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".output", job_name_ + ".cube");

  // padding
  padding_ = options.get(".padding").as<double>();

  // steps
  steps_.y() = options.get(".ysteps").as<Index>();
  steps_.x() = options.get(".xsteps").as<Index>();
  steps_.z() = options.get(".zsteps").as<Index>();

  state_ = options.get(".state").as<QMState>();
  dostateonly_ = options.get(".diff2gs").as<bool>();

  mode_ = options.get(".mode").as<std::string>();
  if (mode_ == "subtract") {
    infile1_ = options.get(".infile1").as<std::string>();
    infile2_ = options.get(".infile2").as<std::string>();
  }
}

void GenCube::calculateCube() {

  XTP_LOG(Log::error, log_)
      << "Reading serialized QM data from " << orbfile_ << std::flush;

  Orbitals orbitals;
  orbitals.ReadFromCpt(orbfile_);

  CubeFile_Writer writer(steps_, padding_, log_);
  XTP_LOG(Log::error, log_) << "Created cube grid" << std::flush;
  writer.WriteFile(output_file_, orbitals, state_, dostateonly_);
  XTP_LOG(Log::error, log_)
      << "Wrote cube data to " << output_file_ << std::flush;
  return;
}

void GenCube::subtractCubes() {

  // open infiles for reading
  std::ifstream in1;
  XTP_LOG(Log::error, log_)
      << " Reading first cube from " << infile1_ << std::flush;
  in1.open(infile1_, std::ios::in);
  std::ifstream in2;
  XTP_LOG(Log::error, log_)
      << " Reading second cube from " << infile2_ << std::flush;
  in2.open(infile2_, std::ios::in);
  std::string s;

  std::ofstream out(output_file_);

  // first two lines of header are garbage
  tools::getline(in1, s);
  out << s << "\n";
  tools::getline(in1, s);
  out << s << " substraction\n";
  tools::getline(in2, s);
  tools::getline(in2, s);

  // read rest from header
  Index natoms;
  double xstart;
  double ystart;
  double zstart;
  // first line
  in1 >> natoms;
  bool do_amplitude = false;
  if (natoms < 0) {
    do_amplitude = true;
  }
  in1 >> xstart;
  in1 >> ystart;
  in1 >> zstart;
  // check from second file
  Index tempint;
  double tempdouble;
  in2 >> tempint;
  if (tempint != natoms) {
    throw std::runtime_error("Atom numbers do not match");
  }
  in2 >> tempdouble;
  if (tempdouble != xstart) {
    throw std::runtime_error("Xstart does not match");
  }
  in2 >> tempdouble;
  if (tempdouble != ystart) {
    throw std::runtime_error("Ystart does not match");
  }
  in2 >> tempdouble;
  if (tempdouble != zstart) {
    throw std::runtime_error("Zstart does not match");
  }

  out << boost::format("%1$lu %2$f %3$f %4$f \n") % natoms % xstart % ystart %
             zstart;

  // grid information from first cube
  double xincr;
  double yincr;
  double zincr;
  Index xsteps;
  Index ysteps;
  Index zsteps;
  in1 >> xsteps;
  in1 >> xincr;
  in1 >> tempdouble;
  in1 >> tempdouble;
  in1 >> ysteps;
  in1 >> tempdouble;
  in1 >> yincr;
  in1 >> tempdouble;
  in1 >> zsteps;
  in1 >> tempdouble;
  in1 >> tempdouble;
  in1 >> zincr;

  // check second cube
  in2 >> tempint;
  if (tempint != xsteps) {
    throw std::runtime_error("xsteps does not match");
  }
  in2 >> tempdouble;
  if (tempdouble != xincr) {
    throw std::runtime_error("xincr does not match");
  }
  in2 >> tempdouble;
  in2 >> tempdouble;
  in2 >> tempint;
  if (tempint != ysteps) {
    throw std::runtime_error("ysteps does not match");
  }
  in2 >> tempdouble;
  in2 >> tempdouble;
  if (tempdouble != yincr) {
    throw std::runtime_error("yincr does not match");
  }
  in2 >> tempdouble;
  in2 >> tempint;
  if (tempint != zsteps) {
    throw std::runtime_error("zsteps does not match");
  }
  in2 >> tempdouble;
  in2 >> tempdouble;
  in2 >> tempdouble;
  if (tempdouble != zincr) {
    throw std::runtime_error("zincr does not match");
  }

  out << boost::format("%1$d %2$f 0.0 0.0 \n") % xsteps % xincr;
  out << boost::format("%1$d 0.0 %2$f 0.0 \n") % ysteps % yincr;
  out << boost::format("%1$d 0.0 0.0 %2$f \n") % zsteps % zincr;

  // atom information

  for (Index iatom = 0; iatom < std::abs(natoms); iatom++) {
    // get center coordinates in Bohr
    double x;
    double y;
    double z;
    Index atnum;
    double crg;

    // get from first cube
    in1 >> atnum;
    in1 >> crg;
    in1 >> x;
    in1 >> y;
    in1 >> z;

    // check second cube
    in2 >> tempint;
    if (tempint != atnum) {
      throw std::runtime_error("atnum does not match");
    }
    in2 >> tempdouble;
    if (tempdouble != crg) {
      throw std::runtime_error("crg does not match");
    }
    in2 >> tempdouble;
    if (tempdouble != x) {
      throw std::runtime_error("x does not match");
    }
    in2 >> tempdouble;
    if (tempdouble != y) {
      throw std::runtime_error("y does not match");
    }
    in2 >> tempdouble;
    if (tempdouble != z) {
      throw std::runtime_error("z does not match");
    }
    out << boost::format("%1$d %2$f %3$f %4$f %5$f\n") % atnum % crg % x % y %
               z;
  }

  if (do_amplitude) {
    Index ntotal;
    Index nis;
    in1 >> ntotal;
    in1 >> nis;

    in2 >> tempint;
    if (tempint != ntotal) {
      throw std::runtime_error("ntotal does not match");
    }
    in2 >> tempint;
    if (tempint != nis) {
      throw std::runtime_error("nis does not match");
    }
    out << boost::format("  1 %1$d \n") % nis;
  }
  // now read data
  double val1;
  double val2;
  for (Index ix = 0; ix < xsteps; ix++) {
    for (Index iy = 0; iy < ysteps; iy++) {
      Index Nrecord = 0;
      for (Index iz = 0; iz < zsteps; iz++) {
        Nrecord++;
        in1 >> val1;
        in2 >> val2;
        if (Nrecord == 6 || iz == zsteps - 1) {
          out << boost::format("%1$E \n") % (val1 - val2);
          Nrecord = 0;
        } else {
          out << boost::format("%1$E ") % (val1 - val2);
        }
      }
    }
  }

  out.close();
  XTP_LOG(Log::error, log_)
      << "Wrote subtracted cube data to " << output_file_ << std::flush;
}

bool GenCube::Run() {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  // calculate new cube
  if (mode_ == "new") {
    calculateCube();
  } else if (mode_ == "subtract") {
    subtractCubes();
  }

  return true;
}
}  // namespace xtp
}  // namespace votca
