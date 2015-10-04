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

#include <fstream>
#include <iostream>
#include <string>

#include "hdf5.h"

namespace votca {  // NOLINT
namespace csg {
using namespace votca::tools;  // NOLINT

/**
    \brief class for reading H5MD trajectory.

    This class implements the H5MD trajectory reading function. The format of the H5MD
    file is defined in
    Pierre de Buyl, Peter H. Colberg, Felix Höfling, H5MD: A structured, efficient, and portable file
    format for molecular data, http://dx.doi.org/10.1016/j.cpc.2014.01.018    The current reference is available here: http://nongnu.org/h5md/
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
  /// Reads dataset that contains vectors.
  template <typename T1>
  T1* ReadVectorData(hid_t ds, hid_t ds_data_type, int row) {
    hsize_t offset[3] = {row, 0, 0};
    hid_t dsp = H5Dget_space(ds);
    H5Sselect_hyperslab(dsp, H5S_SELECT_SET, chunk_rows_, offset, NULL, NULL);
    hid_t mspace1 = H5Screate_simple(vec_components_, chunk_rows_, NULL);
    T1 *data_out = new T1[N_particles_ * vec_components_];
    herr_t status = H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out);
    if (status > 0)
      return data_out;
    else
      throw std::runtime_error("Error");
  }

  /// Reads dataset with scalar values
  template <typename T1>
  T1* ReadScalarData(hid_t ds, hid_t ds_data_type, int row) {
    hsize_t offset[2] = {row, 0};
    hsize_t ch_rows[2] = {1, N_particles_};
    hid_t dsp = H5Dget_space(ds);
    H5Sselect_hyperslab(dsp, H5S_SELECT_SET, ch_rows, offset, NULL, NULL);
    hid_t mspace1 = H5Screate_simple(2, ch_rows, NULL);
    T1 *data_out = new T1[N_particles_];
    herr_t status = H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out);
    if (status > 0)
      return data_out;
    else
      throw std::runtime_error("Error");
  }

  hid_t file_id_;
  hid_t ds_atom_position_;
  hid_t ds_atom_force_;
  hid_t ds_atom_velocity_;
  hid_t ds_atom_id_;

  hid_t particle_group_;
  hid_t atom_position_group_;
  hid_t atom_force_group_;
  hid_t atom_velocity_group_;
  hid_t atom_id_group_;

  int rank_;

  string fname_;
  bool first_frame_;
  bool has_velocity_;
  bool has_force_;
  bool has_id_group_;

  int idx_frame_;

  double *time_set_;
  int *step_set_;

  int N_particles_;
  int vec_components_;
  hsize_t chunk_rows_[3];
  void CheckError(hid_t hid, std::string error_message) {
    if (hid < 0) {
      std::cout << error_message << std::endl;
      throw std::runtime_error(error_message);
    }
  }
};

}  // namespace csg
}  // namespace votca  // NOLINT

#endif  // SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
