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

#include "h5mdtrajectoryreader.h"
#include "hdf5.h"

#include <vector>

namespace votca {
namespace csg {

H5MDTrajectoryReader::H5MDTrajectoryReader(){ 
  has_velocity_ = false;
  has_force_ = false;
  has_id_group_ = false;
}

H5MDTrajectoryReader::~H5MDTrajectoryReader() {
  H5Fclose(file_id_);
}

bool H5MDTrajectoryReader::Open(const string &file) {
  // Checks if we deal with hdf5 file.
  file_id_ = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Check the version of the file.
  hid_t g_h5md = H5Gopen(file_id_, "h5md", H5P_DEFAULT);
  CheckError(g_h5md, "Unable to open /h5md group.");

  hid_t at_version = H5Aopen(g_h5md, "version", H5P_DEFAULT);
  CheckError(at_version, "Unable to read version attribute.");
  int version[2];
  H5Aread(at_version, H5Aget_type(at_version), &version);
  if (version[0] != 1 || (version[0] == 1 && version[1] > 0)) {
    H5Fclose(file_id_);
    std::cout << "Major version " << version[0] << std::endl;
    std::cout << "Minor version " << version[1] << std::endl;
    throw std::ios_base::failure("Wrong version of H5MD file.");
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
    H5Fclose(file_id_);
    throw std::ios_base::failure("The particles group is empty.");
  }

  has_velocity_ = false;
  has_force_ = false;
  first_frame_ = true;

  return true;
}

void H5MDTrajectoryReader::Close() {
  H5Fclose(file_id_);
}

bool H5MDTrajectoryReader::FirstFrame(Topology &top) {  // NOLINT const reference
  herr_t error_id;
  // First frame has access to topology so here we will check if the particle_group exists.
  if (first_frame_) {
    first_frame_ = false;

    std::string *particle_group_name = new std::string(top.getParticleGroup());
    if(*particle_group_name == "")
      throw std::ios_base::failure(
        "Missing particle group in topology. Please set h5md_particle_group tag with name \
        attribute set to the particle group.");
    std::string gr_name = *particle_group_name + "/position";
    atom_position_group_ = H5Gopen(particle_group_, gr_name.c_str(), H5P_DEFAULT);
    CheckError(atom_position_group_, "Unable to open " + gr_name + " group");

    idx_frame_ = -1;
    ds_atom_position_ = H5Dopen(atom_position_group_, "value", H5P_DEFAULT);
    CheckError(ds_atom_position_, "Unable to open " + gr_name + "/value dataset");

    // Reads the box information.
    gr_name = *particle_group_name + "/box";
    hid_t g_box = H5Gopen(particle_group_, gr_name.c_str(), H5P_DEFAULT);
    CheckError(g_box, "Unable to open " + gr_name + " group");
    hid_t at_box_dimension = H5Aopen(g_box, "dimension", H5P_DEFAULT);
    CheckError(at_box_dimension, "Unable to open dimension attribute.");
    int dimension;
    H5Aread(at_box_dimension, H5Aget_type(at_box_dimension), &dimension);
    if (dimension != 3 ) {
      throw std::ios_base::failure("Wrong dimension " + boost::lexical_cast<string>(dimension));
    }

    // Gets the force group.
    std::string force_group_name = *particle_group_name + "/force";
    atom_force_group_ = H5Gopen(particle_group_, force_group_name.c_str(), H5P_DEFAULT);
    if (atom_force_group_ > 0) {
      ds_atom_force_ = H5Dopen(atom_force_group_, "value", H5P_DEFAULT);
      has_force_ = true;
      std::cout << "H5MD: has /force" << std::endl;
    } else {
      has_force_ = false;
    }
    // Gets the velocity group.
    std::string velocity_group_name = *particle_group_name + "/velocity";
    atom_velocity_group_ = H5Gopen(particle_group_, velocity_group_name.c_str(), H5P_DEFAULT);
    if (atom_velocity_group_ > 0) {
      ds_atom_velocity_ = H5Dopen(atom_velocity_group_, "value", H5P_DEFAULT);
      has_velocity_ = true;
      std::cout << "H5MD: has /velocity" << std::endl;
    } else {
      has_velocity_ = false;
    }
    // Gets the id group so that the atom id is taken from this group.
    std::string id_group_name = *particle_group_name  + "/id";
    atom_id_group_ = H5Gopen(particle_group_, id_group_name.c_str(), H5P_DEFAULT);
    if (atom_id_group_ > 0) {
      ds_atom_id_ = H5Dopen(atom_id_group_, "value", H5P_DEFAULT);
      has_id_group_ = true;
      std::cout << "H5MD: has /id group" << std::endl;
    } else {
      has_id_group_ = false;
    }

    // Gets number of particles and dimensions.
    hid_t fs_atom_position_ = H5Dget_space(ds_atom_position_);
    CheckError(fs_atom_position_, "Unable to open atom position space.");
    hsize_t dims[3];
    rank_ = H5Sget_simple_extent_dims(fs_atom_position_, dims, NULL);
    N_particles_ = dims[1];
    vec_components_ = dims[2];

    if (!has_id_group_ && N_particles_ != top.BeadCount()) {
      std::cout << "Warning: The number of beads (" << N_particles_ << ")";
      std::cout << " in the trajectory is different than defined in the topology (" << top.BeadCount() << ")" << std::endl;
      std::cout << "The number of beads from topology will be used!" << std::endl;
      N_particles_ = top.BeadCount();
    }
    chunk_rows_[0] = 1;
    chunk_rows_[1] = N_particles_;
    chunk_rows_[2] = vec_components_;

    delete particle_group_name;
  }
  NextFrame(top);
  return true;
}

/// Reading the data.
bool H5MDTrajectoryReader::NextFrame(Topology &top) {  // NOLINT const reference
  // Reads the position row.
  idx_frame_++;
  std::cout << '\r' << "Reading frame: " << idx_frame_;
  std::cout.flush();
  double *positions;
  try {
    positions = ReadVectorData<double>(ds_atom_position_, H5T_NATIVE_DOUBLE, idx_frame_);
  } catch (const std::runtime_error &e) {
    return false;
  }

  double *forces = NULL;
  double *velocities = NULL;
  int *ids = NULL;
  if (has_velocity_) {
    velocities = ReadVectorData<double>(ds_atom_velocity_, H5T_NATIVE_DOUBLE, idx_frame_);
  }

  if (has_force_) {
    forces = ReadVectorData<double>(ds_atom_force_, H5T_NATIVE_DOUBLE, idx_frame_);
  }

  if (has_id_group_) {
    ids = ReadScalarData<int>(ds_atom_id_, H5T_NATIVE_INT, idx_frame_);
  }

  //Process atoms.
  for (int at_idx = 0; at_idx < N_particles_; at_idx++) {
    double x, y, z;
    int array_index = at_idx*vec_components_;
    x = positions[array_index];
    y = positions[array_index + 1];
    z = positions[array_index + 2];
    // Set atom id, or it is an index of a row in dataset or from id dataset.
    int atom_id = at_idx;
    if (has_id_group_) {
      if (ids[at_idx] == -1)  // ignore values where id == -1
        continue;
      atom_id = ids[at_idx] - 1;
    }

    // Topology has to be defined in the xml file or in other
    // topology files. The h5md only stores the trajectory data.
    Bead *b = top.getBead(atom_id);
    b->setPos(vec(x, y, z));
    if (has_velocity_) {
      double vx, vy, vz;
      vx = velocities[array_index];
      vy = velocities[array_index + 1];
      vz = velocities[array_index + 2];
      b->setVel(vec(vx, vy, vz));
    }

    if (has_force_) {
      double fx, fy, fz;
      fx = forces[array_index];
      fy = forces[array_index + 1];
      fz = forces[array_index + 2];
      b->setF(vec(fx, fy, fz));
    }
  }

  // Clean up pointers.
  delete[] positions;
  delete[] forces;
  delete[] velocities;
  delete[] ids;

  return true;
}

}  // end csg
}  // end votca

