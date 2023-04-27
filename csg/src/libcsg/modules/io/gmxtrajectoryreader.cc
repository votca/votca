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

// Standard includes
#include <cstdlib>
#include <filesystem>
#include <iostream>

// Third party includes
#include <gromacs/utility/programcontext.h>
#include <gromacs/version.h>

// Local VOTCA includes
#include "votca/csg/topology.h"

// Local private VOTCA includes
#include "gmxtrajectoryreader.h"

namespace votca {
namespace csg {

bool GMXTrajectoryReader::Open(const std::string &file) {
  filename_ = file;
  return true;
}

void GMXTrajectoryReader::Close() { close_trx(gmx_status_); }

bool GMXTrajectoryReader::FirstFrame(Topology &conf) {
  gmx_output_env_t *oenv;
#if GMX_VERSION >= 20210000
  gmx::TimeUnit timeunit = gmx::TimeUnit::Picoseconds;
  XvgFormat xvg = XvgFormat::None;
#else
  time_unit_t timeunit = time_ps;
  xvg_format_t xvg = exvgNONE;
#endif
#if GMX_VERSION >= 20230000
  const std::filesystem::path path = filename_;
#else
  char *path = (char *)filename_.c_str();
#endif
  output_env_init(&oenv, gmx::getProgramContext(), timeunit, FALSE, xvg, 0);
  if (!read_first_frame(oenv, &gmx_status_, path, &gmx_frame_,
                        TRX_READ_X | TRX_READ_V | TRX_READ_F)) {
    throw std::runtime_error(std::string("cannot open ") + filename_);
  }
  output_env_done(oenv);

  Eigen::Matrix3d m;
  for (Index i = 0; i < 3; i++) {
    for (Index j = 0; j < 3; j++) {
      m(i, j) = gmx_frame_.box[j][i];
    }
  }
  conf.setBox(m);
  conf.setTime(gmx_frame_.time);
  conf.setStep(gmx_frame_.step);
  std::cout << std::endl;

  if (gmx_frame_.natoms != (Index)conf.Beads().size()) {
    throw std::runtime_error(
        "number of beads in trajectory do not match topology");
  }

  // conf.HasPos(true);
  // conf.HasF( gmx_frame_.bF);

  for (Index i = 0; i < gmx_frame_.natoms; i++) {
    Eigen::Vector3d r = {gmx_frame_.x[i][XX], gmx_frame_.x[i][YY],
                         gmx_frame_.x[i][ZZ]};
    conf.getBead(i)->setPos(r);
    if (gmx_frame_.bF) {
      Eigen::Vector3d f = {gmx_frame_.f[i][XX], gmx_frame_.f[i][YY],
                           gmx_frame_.f[i][ZZ]};
      conf.getBead(i)->setF(f);
    }
    if (gmx_frame_.bV) {
      Eigen::Vector3d v = {gmx_frame_.v[i][XX], gmx_frame_.v[i][YY],
                           gmx_frame_.v[i][ZZ]};
      conf.getBead(i)->setVel(v);
    }
  }
  return true;
}

bool GMXTrajectoryReader::NextFrame(Topology &conf) {
  gmx_output_env_t *oenv;
#if GMX_VERSION >= 20210000
  gmx::TimeUnit timeunit = gmx::TimeUnit::Picoseconds;
  XvgFormat xvg = XvgFormat::None;
#else
  time_unit_t timeunit = time_ps;
  xvg_format_t xvg = exvgNONE;
#endif
  output_env_init(&oenv, gmx::getProgramContext(), timeunit, FALSE, xvg, 0);
  if (!read_next_frame(oenv, gmx_status_, &gmx_frame_)) {
    return false;
  }
  output_env_done(oenv);

  Eigen::Matrix3d m;
  for (Index i = 0; i < 3; i++) {
    for (Index j = 0; j < 3; j++) {
      m(i, j) = gmx_frame_.box[j][i];
    }
  }
  conf.setTime(gmx_frame_.time);
  conf.setStep(gmx_frame_.step);
  conf.setBox(m);

  // conf.HasF( gmx_frame_.bF);

  for (Index i = 0; i < gmx_frame_.natoms; i++) {
    Eigen::Vector3d r = {gmx_frame_.x[i][XX], gmx_frame_.x[i][YY],
                         gmx_frame_.x[i][ZZ]};
    conf.getBead(i)->setPos(r);
    if (gmx_frame_.bF) {
      Eigen::Vector3d f = {gmx_frame_.f[i][XX], gmx_frame_.f[i][YY],
                           gmx_frame_.f[i][ZZ]};
      conf.getBead(i)->setF(f);
    }
    if (gmx_frame_.bV) {
      Eigen::Vector3d v = {gmx_frame_.v[i][XX], gmx_frame_.v[i][YY],
                           gmx_frame_.v[i][ZZ]};
      conf.getBead(i)->setVel(v);
    }
  }
  return true;
}

}  // namespace csg
}  // namespace votca
