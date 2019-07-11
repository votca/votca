#include <memory>

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

#include "h5mdtrajectoryreader.h"
#include "hdf5.h"

#include <vector>

namespace votca {
namespace csg {

using namespace std;

H5MDTrajectoryReader::H5MDTrajectoryReader() {
  has_velocity_ = H5MDTrajectoryReader::NONE;
  has_force_ = H5MDTrajectoryReader::NONE;
  has_id_group_ = H5MDTrajectoryReader::NONE;
  has_box_ = H5MDTrajectoryReader::NONE;
}

H5MDTrajectoryReader::~H5MDTrajectoryReader() {
  if (file_opened_) {
    H5Fclose(file_id_);
    file_opened_ = false;
  }
}

bool H5MDTrajectoryReader::Open(const string &file) {
  // Checks if we deal with hdf5 file.
  if (!H5Fis_hdf5(file.c_str())) {
    cout << file << " is not recognise as HDF5 file format" << endl;
    return false;
  }
  file_id_ = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  file_opened_ = true;

  // Check the version of the file.
  hid_t g_h5md = H5Gopen(file_id_, "h5md", H5P_DEFAULT);
  CheckError(g_h5md, "Unable to open /h5md group.");

  hid_t at_version = H5Aopen(g_h5md, "version", H5P_DEFAULT);
  CheckError(at_version, "Unable to read version attribute.");
  int version[2];
  H5Aread(at_version, H5Aget_type(at_version), &version);
  if (version[0] != 1 || (version[0] == 1 && version[1] > 1)) {
    cout << "Found H5MD version: " << version[0] << "." << version[1] << endl;
    throw ios_base::failure("Wrong version of H5MD file.");
  }

  // Clean up.
  H5Aclose(at_version);
  H5Gclose(g_h5md);

  // Checks if the particles group exists and what is the number of members.
  particle_group_ = H5Gopen(file_id_, "particles", H5P_DEFAULT);
  CheckError(particle_group_, "Unable to open /particles group.");
  hsize_t num_obj = 0;
  H5Gget_num_objs(particle_group_, &num_obj);
  if (num_obj == 0) {
    throw ios_base::failure("The particles group is empty.");
  }

  first_frame_ = true;

  // Handle errors by internal check up.
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  return true;
}

void H5MDTrajectoryReader::Close() {
  if (file_opened_) {
    H5Fclose(file_id_);
    file_opened_ = false;
  }
}

void H5MDTrajectoryReader::Initialize(Topology &top) {
  string particle_group_name_ = top.getParticleGroup();
  if (particle_group_name_.compare("unassigned") == 0)
    throw ios_base::failure(
        "Missing particle group in topology. Please set `h5md_particle_group` "
        "tag with `name` attribute set to the particle group.");
  string position_group_name = particle_group_name_ + "/position";
  atom_position_group_ =
      H5Gopen(particle_group_, position_group_name.c_str(), H5P_DEFAULT);
  CheckError(atom_position_group_,
             "Unable to open " + position_group_name + " group");

  idx_frame_ = -1;
  ds_atom_position_ = H5Dopen(atom_position_group_, "value", H5P_DEFAULT);
  CheckError(ds_atom_position_,
             "Unable to open " + position_group_name + "/value dataset");

  // Reads the box information.
  string box_gr_name = particle_group_name_ + "/box";
  hid_t g_box = H5Gopen(particle_group_, box_gr_name.c_str(), H5P_DEFAULT);
  CheckError(g_box, "Unable to open " + box_gr_name + " group");
  hid_t at_box_dimension = H5Aopen(g_box, "dimension", H5P_DEFAULT);
  CheckError(at_box_dimension, "Unable to open dimension attribute.");
  int dimension;
  H5Aread(at_box_dimension, H5Aget_type(at_box_dimension), &dimension);
  if (dimension != 3) {
    throw ios_base::failure("Wrong dimension " +
                            boost::lexical_cast<string>(dimension));
  }
  // TODO: check if boundary is periodic.
  string box_edges_name = particle_group_name_ + "/box/edges";
  if (GroupExists(particle_group_, box_edges_name)) {
    g_box = H5Gopen(particle_group_, box_gr_name.c_str(), H5P_DEFAULT);
    edges_group_ = H5Gopen(g_box, "edges", H5P_DEFAULT);
    ds_edges_group_ = H5Dopen(edges_group_, "value", H5P_DEFAULT);
    cout << "H5MD: has /box/edges" << endl;
    cout << "H5MD: time dependent box size" << endl;
    has_box_ = H5MDTrajectoryReader::TIMEDEPENDENT;
  } else {
    cout << "H5MD: static box" << endl;
    hid_t ds_edges = H5Dopen(g_box, "edges", H5P_DEFAULT);
    CheckError(ds_edges, "Unable to open /box/edges");
    unique_ptr<double[]> box = unique_ptr<double[]>{new double[3]};
    ReadStaticData<double[]>(ds_edges, H5T_NATIVE_DOUBLE, box);
    cout << "H5MD: Found box " << box[0] << " x " << box[1] << " x " << box[2]
         << endl;
    // Sets box size.
    m = Eigen::Matrix3d::Zero();
    m(0, 0) = box[0];
    m(1, 1) = box[1];
    m(2, 2) = box[2];
    top.setBox(m);
    has_box_ = H5MDTrajectoryReader::STATIC;
  }
  H5Gclose(g_box);

  // Gets the force group.
  string force_group_name = particle_group_name_ + "/force";
  if (GroupExists(particle_group_, force_group_name)) {
    atom_force_group_ =
        H5Gopen(particle_group_, force_group_name.c_str(), H5P_DEFAULT);
    ds_atom_force_ = H5Dopen(atom_force_group_, "value", H5P_DEFAULT);
    has_force_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    cout << "H5MD: has /force" << endl;
  } else {
    has_force_ = H5MDTrajectoryReader::NONE;
  }

  // Gets the velocity group.
  string velocity_group_name = particle_group_name_ + "/velocity";
  if (GroupExists(particle_group_, velocity_group_name)) {
    atom_velocity_group_ =
        H5Gopen(particle_group_, velocity_group_name.c_str(), H5P_DEFAULT);
    ds_atom_velocity_ = H5Dopen(atom_velocity_group_, "value", H5P_DEFAULT);
    has_velocity_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    cout << "H5MD: has /velocity" << endl;
  } else {
    has_velocity_ = H5MDTrajectoryReader::NONE;
  }

  // Gets the id group so that the atom id is taken from this group.
  string id_group_name = particle_group_name_ + "/id";
  if (GroupExists(particle_group_, id_group_name)) {
    atom_id_group_ =
        H5Gopen(particle_group_, id_group_name.c_str(), H5P_DEFAULT);
    ds_atom_id_ = H5Dopen(atom_id_group_, "value", H5P_DEFAULT);
    has_id_group_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    cout << "H5MD: has /id group" << endl;
  } else {
    has_id_group_ = H5MDTrajectoryReader::NONE;
  }

  // Gets number of particles and dimensions.
  hid_t fs_atom_position_ = H5Dget_space(ds_atom_position_);
  CheckError(fs_atom_position_, "Unable to open atom position space.");
  hsize_t dims[3];
  rank_ = H5Sget_simple_extent_dims(fs_atom_position_, dims, NULL);
  N_particles_ = dims[1];
  vec_components_ = dims[2];
  max_idx_frame_ = dims[0] - 1;

  // TODO: reads mass, charge and particle type.

  if (has_id_group_ == H5MDTrajectoryReader::NONE && top.BeadCount() > 0 &&
      N_particles_ != top.BeadCount()) {
    cout << "Warning: The number of beads (" << N_particles_ << ")";
    cout << " in the trajectory is different than defined in the topology ("
         << top.BeadCount() << ")" << endl;
    cout << "The number of beads from topology will be used!" << endl;
    N_particles_ = top.BeadCount();
  }
}

bool H5MDTrajectoryReader::FirstFrame(Topology &top) {  // NOLINT const
  // reference
  if (first_frame_) {
    first_frame_ = false;
    Initialize(top);
  }
  NextFrame(top);
  return true;
}

/// Reading the data.
bool H5MDTrajectoryReader::NextFrame(Topology &top) {  // NOLINT const reference
  // Reads the position row.
  idx_frame_++;
  if (idx_frame_ > max_idx_frame_) return false;

  cout << '\r' << "Reading frame: " << idx_frame_ << "\n";
  cout.flush();
  // Set volume of box because top on workers somehow does not have this
  // information.
  if (has_box_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
    unique_ptr<double[]> box = unique_ptr<double[]>{new double[3]};
    ReadBox(ds_edges_group_, H5T_NATIVE_DOUBLE, idx_frame_, box);
    m = Eigen::Matrix3d::Zero();
    m(0, 0) = box.get()[0];
    m(1, 1) = box.get()[1];
    m(2, 2) = box.get()[2];
    cout << "Time dependent box:" << endl;
    cout << m << endl;
  }
  top.setBox(m);

  double *positions;
  double *forces = NULL;
  double *velocities = NULL;
  int *ids = NULL;

  try {
    positions = ReadVectorData<double>(ds_atom_position_, H5T_NATIVE_DOUBLE,
                                       idx_frame_);
  } catch (const runtime_error &e) {
    return false;
  }

  if (has_velocity_ != H5MDTrajectoryReader::NONE) {
    velocities = ReadVectorData<double>(ds_atom_velocity_, H5T_NATIVE_DOUBLE,
                                        idx_frame_);
  }

  if (has_force_ != H5MDTrajectoryReader::NONE) {
    forces =
        ReadVectorData<double>(ds_atom_force_, H5T_NATIVE_DOUBLE, idx_frame_);
  }

  if (has_id_group_ != H5MDTrajectoryReader::NONE) {
    ids = ReadScalarData<int>(ds_atom_id_, H5T_NATIVE_INT, idx_frame_);
  }

  // Process atoms.
  for (int at_idx = 0; at_idx < N_particles_; at_idx++) {
    double x, y, z;
    int array_index = at_idx * vec_components_;
    x = positions[array_index];
    y = positions[array_index + 1];
    z = positions[array_index + 2];
    // Set atom id, or it is an index of a row in dataset or from id dataset.
    int atom_id = at_idx;
    if (has_id_group_ != H5MDTrajectoryReader::NONE) {
      if (ids[at_idx] == -1)  // ignore values where id == -1
        continue;
      atom_id = ids[at_idx] - 1;
    }

    // Topology has to be defined in the xml file or in other
    // topology files. The h5md only stores the trajectory data.
    Bead *b = top.getBead(atom_id);
    if (b == nullptr)
      throw runtime_error("Bead not found: " +
                          boost::lexical_cast<string>(atom_id));

    b->setPos(Eigen::Vector3d(x, y, z));
    if (has_velocity_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
      double vx, vy, vz;
      vx = velocities[array_index];
      vy = velocities[array_index + 1];
      vz = velocities[array_index + 2];
      b->setVel(Eigen::Vector3d(vx, vy, vz));
    }

    if (has_force_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
      double fx, fy, fz;
      fx = forces[array_index];
      fy = forces[array_index + 1];
      fz = forces[array_index + 2];
      b->setF(Eigen::Vector3d(fx, fy, fz));
    }
  }

  // Clean up pointers.
  delete[] positions;
  if (has_force_ == H5MDTrajectoryReader::TIMEDEPENDENT) delete[] forces;
  if (has_velocity_ == H5MDTrajectoryReader::TIMEDEPENDENT) delete[] velocities;
  if (has_id_group_ == H5MDTrajectoryReader::TIMEDEPENDENT) delete[] ids;

  return true;
}

void H5MDTrajectoryReader::ReadBox(hid_t ds, hid_t ds_data_type, int row,
                                   unique_ptr<double[]> &data_out) {
  hsize_t offset[2];
  offset[0] = row;
  offset[1] = 0;
  hsize_t ch_rows[2];
  ch_rows[0] = 1;
  ch_rows[1] = 3;
  hid_t dsp = H5Dget_space(ds);
  H5Sselect_hyperslab(dsp, H5S_SELECT_SET, offset, NULL, ch_rows, NULL);
  hid_t mspace1 = H5Screate_simple(2, ch_rows, NULL);
  herr_t status =
      H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out.get());
  if (status < 0) {
    throw std::runtime_error("Error ReadScalarData: " +
                             boost::lexical_cast<std::string>(status));
  }
}

}  // namespace csg
}  // namespace votca
