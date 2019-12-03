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

#ifndef SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
#define SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_

#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "hdf5.h"

namespace votca {  // NOLINT
namespace csg {

/**
    \brief class for reading H5MD trajectory.

    This class implements the H5MD trajectory reading function. The format of
   the H5MD file is defined in Pierre de Buyl, Peter H. Colberg, Felix Höfling,
   H5MD: A structured, efficient, and portable file format for molecular data,
   http://dx.doi.org/10.1016/j.cpc.2014.01.018 The current reference is
   available here: http://nongnu.org/h5md/
*/
class H5MDTrajectoryReader : public TrajectoryReader {
 public:
  H5MDTrajectoryReader();
  ~H5MDTrajectoryReader() override;

  /// Opens original trajectory file.
  bool Open(const std::string &file) override;

  /// Initialize data structures.
  void Initialize(Topology &top);

  /// Reads in the first frame.
  bool FirstFrame(Topology &conf) override;

  /// Reads in the next frame.
  bool NextFrame(Topology &conf) override;

  /// Closes original trajectory file.
  void Close() override;

 private:
  enum DatasetState { NONE, STATIC, TIMEDEPENDENT };

  /// Reads dataset that contains vectors.
  template <typename T1>
  T1 *ReadVectorData(hid_t ds, hid_t ds_data_type, Index row) {
    hsize_t offset[3];
    offset[0] = hsize_t(row);
    offset[1] = 0;
    offset[2] = 0;
    hsize_t chunk_rows[3];
    chunk_rows[0] = 1;
    chunk_rows[1] = N_particles_;
    chunk_rows[2] = vec_components_;
    hid_t dsp = H5Dget_space(ds);
    H5Sselect_hyperslab(dsp, H5S_SELECT_SET, offset, nullptr, chunk_rows,
                        nullptr);
    hid_t mspace1 = H5Screate_simple(vec_components_, chunk_rows, nullptr);
    T1 *data_out = new T1[N_particles_ * vec_components_];
    herr_t status =
        H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out);
    if (status < 0) {
      throw std::runtime_error("Error ReadVectorData: " +
                               boost::lexical_cast<std::string>(status));
    } else {
      return data_out;
    }
  }

  /// Reads dataset with scalar values.
  template <typename T1>
  T1 *ReadScalarData(hid_t ds, hid_t ds_data_type, Index row) {
    hsize_t offset[2];
    offset[0] = row;
    offset[1] = 0;
    hsize_t ch_rows[2];
    ch_rows[0] = 1;
    ch_rows[1] = N_particles_;
    hid_t dsp = H5Dget_space(ds);
    H5Sselect_hyperslab(dsp, H5S_SELECT_SET, offset, nullptr, ch_rows, nullptr);
    hid_t mspace1 = H5Screate_simple(2, ch_rows, nullptr);
    T1 *data_out = new T1[N_particles_];
    herr_t status =
        H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out);
    if (status < 0) {
      throw std::runtime_error("Error ReadScalarData: " +
                               boost::lexical_cast<std::string>(status));
    } else {
      return data_out;
    }
  }

  template <typename T1>
  void ReadStaticData(hid_t ds, hid_t ds_data_type,
                      std::unique_ptr<T1> &outbuf) {
    herr_t status =
        H5Dread(ds, ds_data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, outbuf.get());
    if (status < 0) {
      H5Eprint(H5E_DEFAULT, stderr);
    }
  }

  void ReadBox(hid_t ds, hid_t ds_data_type, Index row,
               std::unique_ptr<double[]> &data_out);

  void CheckError(hid_t hid, std::string error_message) {
    if (hid < 0) {
      // H5Eprint(H5E_DEFAULT, stderr);
      throw std::runtime_error(error_message);
    }
  }

  bool GroupExists(hid_t file_id, std::string path) {
    H5G_stat_t info;
    herr_t status = H5Gget_objinfo(file_id, path.c_str(), 0, &info);
    if (status < 0) {
      return false;
    }
    return info.type == H5G_GROUP;
  }

  hid_t file_id_;
  hid_t ds_atom_position_;
  hid_t ds_atom_force_;
  hid_t ds_atom_velocity_;
  hid_t ds_atom_id_;
  hid_t ds_edges_group_;

  hid_t particle_group_;
  hid_t atom_position_group_;
  hid_t atom_force_group_;
  hid_t atom_velocity_group_;
  hid_t atom_id_group_;
  hid_t edges_group_;

  Index rank_;

  std::string fname_;
  bool first_frame_;

  // Flags about datasets.
  DatasetState has_velocity_;
  DatasetState has_force_;
  DatasetState has_id_group_;
  DatasetState has_box_;

  bool file_opened_;

  // Current frame indicator.
  Index idx_frame_;
  Index max_idx_frame_;

  // Number of particles. This is static among time.
  Index N_particles_;
  //
  int vec_components_;

  // Box matrix.
  Eigen::Matrix3d m;
};

}  // namespace csg
}  // namespace votca

#endif  // SRC_LIBCSG_MODULES_IO_H5MDTRAJECTORYREADER_H_
