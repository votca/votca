/*
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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

#ifndef SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
#define SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_

#include <votca/csg/trajectoryreader.h>
#include "H5Cpp.h"

#include <string>
#include <iostream>
#include <fstream>

namespace votca {  // NOLINT
namespace csg {
using namespace votca::tools;  // NOLINT

/**
    \brief class for reading H5MD trajectory.

    This class implements the H5MD trajectory reading function. The format of the H5MD
    file is defined in
    Pierre de Buyl, Peter H. Colberg, Felix Höfling, H5MD: A structured, efficient, and portable file
    format for molecular data, http://dx.doi.org/10.1016/j.cpc.2014.01.018.

    The current reference is available here: http://nongnu.org/h5md/
*/
class H5MDTrajectoryReader : public TrajectoryReader {
 public:
  H5MDTrajectoryReader();
  ~H5MDTrajectoryReader();

  /// Opens original trajectory file.
  bool Open(const string &file);

  /// Reads in the first frame.
  bool FirstFrame(Topology &conf);  // NOLINT

  /// Reads in the next frame.
  bool NextFrame(Topology &conf);  // NOLINT

  /// close original trajectory file.
  void Close();

 private:
  double* ReadData(H5::DataSet *ds, int row);

  H5::H5File *h5file_;
  H5::Group *particle_group_;
  H5::Group *atom_position_group_;
  H5::Group *atom_force_group_;
  H5::Group *atom_velocity_group_;

  H5::DataSet *ds_atom_position_;
  H5::DataSet *ds_atom_force_;
  H5::DataSet *ds_atom_velocity_;

  H5::DataSet *ds_time_;
  H5::DataSet *ds_step_;

  int rank_;

  H5::DataSet *ds_edges_;

  string fname_;
  bool first_frame_;
  bool has_velocity_;
  bool has_force_;
  bool has_static_box_;
  bool has_static_species_;

  int idx_frame_;

  double *time_set_;
  int *step_set_;

  int N_particles_;
  int variables_;
  hsize_t chunk_rows_[3];
};

}  // namespace csg
}  // namespace votca  // NOLINT

#endif  // SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
