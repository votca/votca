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
#include <memory>
#include <string>
#include <vector>

// Third party includes
#include <hdf5.h>
#include <boost/regex.hpp>

// Local private VOTCA includes
#include "h5mdtrajectoryreader.h"

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

bool H5MDTrajectoryReader::Open(const std::string &file) {
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
  int version[2] = {0, 0};
  H5Aread(at_version, H5Aget_type(at_version), &version);
  if (version[0] != 1 || version[1] > 1) {
    cout << "Found H5MD version: " << version[0] << "." << version[1] << endl;
    throw ios_base::failure("Wrong version of H5MD file.");
  }

  // Checks if the particles group exists and what is the number of members.
  particle_group_ = H5Gopen(file_id_, "particles", H5P_DEFAULT);
  CheckError(particle_group_, "Unable to open /particles group.");
  hsize_t num_obj = 0;
  H5Gget_num_objs(particle_group_, &num_obj);
  if (num_obj == 0) {
    throw ios_base::failure("The particles group is empty.");
  }

  // Check if the unit module is enabled.
  hid_t modules = H5Gopen(g_h5md, "modules", H5P_DEFAULT);
  std::cout << "open modules group" << modules << std::endl;
  if (modules > 0) {
    hid_t module_units = H5Gopen(modules, "units", H5P_DEFAULT);
    if (module_units > 0) {
      std::cout << "found units module - position and forces will be scaled " << std::endl;
      unit_module_enabled_ = true;
      H5Gclose(module_units);
    }
    H5Gclose(modules);
  }

  first_frame_ = true;

  // Handle errors by internal check up.
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

  // Clean up.
  H5Aclose(at_version);
  H5Gclose(g_h5md);

  return true;
}

void H5MDTrajectoryReader::Close() {
  if (file_opened_) {
    H5Fclose(file_id_);
    file_opened_ = false;
  }
}

void H5MDTrajectoryReader::Initialize(Topology &top) {
  std::string particle_group_name_ = top.getParticleGroup();
  if (particle_group_name_.compare("unassigned") == 0) {
    throw ios_base::failure(
        "Missing particle group in topology. Please set `h5md_particle_group` "
        "tag with `name` attribute set to the particle group.");
  }
  std::string position_group_name = particle_group_name_ + "/position";
  atom_position_group_ =
      H5Gopen(particle_group_, position_group_name.c_str(), H5P_DEFAULT);
  CheckError(atom_position_group_,
             "Unable to open " + position_group_name + " group");
  idx_frame_ = -1;
  ds_atom_position_ = H5Dopen(atom_position_group_, "value", H5P_DEFAULT);
  CheckError(ds_atom_position_,
             "Unable to open " + position_group_name + "/value dataset");

  // Reads the box information.
  std::string box_gr_name = particle_group_name_ + "/box";
  hid_t g_box = H5Gopen(particle_group_, box_gr_name.c_str(), H5P_DEFAULT);
  CheckError(g_box, "Unable to open " + box_gr_name + " group");
  hid_t at_box_dimension = H5Aopen(g_box, "dimension", H5P_DEFAULT);
  CheckError(at_box_dimension, "Unable to open dimension attribute.");
  int dimension;
  H5Aread(at_box_dimension, H5Aget_type(at_box_dimension), &dimension);
  if (dimension != 3) {
    throw ios_base::failure("Wrong dimension " +
        boost::lexical_cast<std::string>(dimension));
  }
  // TODO: check if boundary is periodic.
  std::string box_edges_name = particle_group_name_ + "/box/edges";
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
    std::unique_ptr<double[]> box = std::unique_ptr<double[]>{new double[3]};
    ReadStaticData<double[]>(ds_edges, H5T_NATIVE_DOUBLE, box);
    cout << "H5MD: Found box " << box[0] << " x " << box[1] << " x " << box[2]
         << endl;
    // Sets box size.
    m = Eigen::Matrix3d::Zero();
    m(0, 0) = box[0] * length_scaling_;
    m(1, 1) = box[1] * length_scaling_;
    m(2, 2) = box[2] * length_scaling_;
    top.setBox(m);
    has_box_ = H5MDTrajectoryReader::STATIC;
  }
  H5Gclose(g_box);

  // Gets the force group.
  std::string force_group_name = particle_group_name_ + "/force";
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
  std::string velocity_group_name = particle_group_name_ + "/velocity";
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
  std::string id_group_name = particle_group_name_ + "/id";
  if (GroupExists(particle_group_, id_group_name)) {
    atom_id_group_ =
        H5Gopen(particle_group_, id_group_name.c_str(), H5P_DEFAULT);
    ds_atom_id_ = H5Dopen(atom_id_group_, "value", H5P_DEFAULT);
    has_id_group_ = H5MDTrajectoryReader::TIMEDEPENDENT;
    cout << "H5MD: has /id group" << endl;
  } else {
    has_id_group_ = H5MDTrajectoryReader::NONE;
  }

  // Reads unit system - if enabled.
  if (unit_module_enabled_) {
    length_scaling_ = ReadScaleFactor(ds_atom_position_, "position");
    force_scaling_ = ReadScaleFactor(ds_atom_force_, "force");
    velocity_scaling_ = ReadScaleFactor(ds_atom_velocity_, "velocity");
  }

  // Gets number of particles and dimensions.
  hid_t fs_atom_position_ = H5Dget_space(ds_atom_position_);
  CheckError(fs_atom_position_, "Unable to open atom position space.");
  hsize_t dims[3];
  rank_ = H5Sget_simple_extent_dims(fs_atom_position_, dims, nullptr);
  N_particles_ = dims[1];
  vec_components_ = (int) dims[2];
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
  if (idx_frame_ > max_idx_frame_) {
    return false;
  }

  cout << '\r' << "Reading frame: " << idx_frame_ << "\n";
  cout.flush();
  // Set volume of box because top on workers somehow does not have this
  // information.
  if (has_box_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
    std::unique_ptr<double[]> box = std::unique_ptr<double[]>{new double[3]};
    ReadBox(ds_edges_group_, H5T_NATIVE_DOUBLE, idx_frame_, box);
    m = Eigen::Matrix3d::Zero();
    m(0, 0) = box.get()[0] * length_scaling_;
    m(1, 1) = box.get()[1] * length_scaling_;
    m(2, 2) = box.get()[2] * length_scaling_;
    cout << "Time dependent box:" << endl;
    cout << m << endl;
  }
  top.setBox(m);

  double *positions;
  double *forces = nullptr;
  double *velocities = nullptr;
  int *ids = nullptr;

  try {
    positions = ReadVectorData<double>(ds_atom_position_, H5T_NATIVE_DOUBLE,
                                       idx_frame_);
  } catch (const std::runtime_error &e) {
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
  for (Index at_idx = 0; at_idx < N_particles_; at_idx++) {
    double x, y, z;
    Index array_index = at_idx * vec_components_;
    x = positions[array_index] * length_scaling_;
    y = positions[array_index + 1] * length_scaling_;
    z = positions[array_index + 2] * length_scaling_;
    // Set atom id, or it is an index of a row in dataset or from id dataset.
    Index atom_id = at_idx;
    if (has_id_group_ != H5MDTrajectoryReader::NONE) {
      if (ids[at_idx] == -1) {  // ignore values where id == -1
        continue;
      }
      atom_id = ids[at_idx];
    }

    // Topology has to be defined in the xml file or in other
    // topology files. The h5md only stores the trajectory data.
    Bead *b = top.getBead(atom_id);
    if (b == nullptr) {
      throw std::runtime_error("Bead not found: " +
          boost::lexical_cast<std::string>(atom_id));
    }

    b->setPos(Eigen::Vector3d(x, y, z));
    if (has_velocity_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
      double vx, vy, vz;
      vx = velocities[array_index] * velocity_scaling_;
      vy = velocities[array_index + 1] * velocity_scaling_;
      vz = velocities[array_index + 2] * velocity_scaling_;
      b->setVel(Eigen::Vector3d(vx, vy, vz));
    }

    if (has_force_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
      double fx, fy, fz;
      fx = forces[array_index] * force_scaling_;
      fy = forces[array_index + 1] * force_scaling_;
      fz = forces[array_index + 2] * force_scaling_;
      b->setF(Eigen::Vector3d(fx, fy, fz));
    }
  }

  // Clean up pointers.
  delete[] positions;
  if (has_force_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
    delete[] forces;
  }
  if (has_velocity_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
    delete[] velocities;
  }
  if (has_id_group_ == H5MDTrajectoryReader::TIMEDEPENDENT) {
    delete[] ids;
  }

  return true;
}

void H5MDTrajectoryReader::ReadBox(hid_t ds, hid_t ds_data_type, Index row,
                                   std::unique_ptr<double[]> &data_out) {
  hsize_t offset[2];
  offset[0] = row;
  offset[1] = 0;
  hsize_t ch_rows[2];
  ch_rows[0] = 1;
  ch_rows[1] = 3;
  hid_t dsp = H5Dget_space(ds);
  H5Sselect_hyperslab(dsp, H5S_SELECT_SET, offset, nullptr, ch_rows, nullptr);
  hid_t mspace1 = H5Screate_simple(2, ch_rows, nullptr);
  herr_t status =
      H5Dread(ds, ds_data_type, mspace1, dsp, H5P_DEFAULT, data_out.get());
  if (status < 0) {
    throw std::runtime_error("Error ReadScalarData: " +
        boost::lexical_cast<std::string>(status));
  }
}

double H5MDTrajectoryReader::ReadScaleFactor(hid_t ds, std::string unit_type) {
  hid_t unit_attr = H5Aopen(ds, "unit", H5P_DEFAULT);
  double scaling_factor = 1.0;
  if (unit_attr > 0) {
    hid_t type_id = H5Aget_type(unit_attr);
    hid_t atype_mem = H5Tget_native_type(type_id, H5T_DIR_ASCEND);
    hsize_t str_size = H5Tget_size(type_id);
    char buffer[80] = {0};  // buffer for attribute
    H5Aread(unit_attr, atype_mem, &buffer);
    H5Tclose(atype_mem);
    H5Tclose(type_id);
    std::string value = std::string(buffer, str_size);

    // Read units
    tools::Tokenizer tok(value, " ");
    for (auto v : tok) {
      boost::smatch suffix_match;
      int unit_pow = 1;
      if (boost::regex_match(v, suffix_match, suffix_units, boost::match_extra)) {
        //If the prefix has numeric suffix then use it in the formula for scaling factor.
        unit_pow = std::stoi(suffix_match[2]);
      }
      auto votca_scaling_factor = length_units_votca_scaling_factors.find(v);
      if (votca_scaling_factor != length_units_votca_scaling_factors.end()) {
        scaling_factor = scaling_factor * pow(votca_scaling_factor->second, unit_pow);
      }
    }
    std::cout << "Found units " << value << " for " << unit_type << ".";
    if (scaling_factor != 1.0) {
      std::cout << " Values scaled by a factor " << scaling_factor;
    }
    std::cout << std::endl;
  }

  return scaling_factor;
}

}  // namespace csg
}  // namespace votca
