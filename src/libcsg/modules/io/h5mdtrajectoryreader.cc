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

H5MDTrajectoryReader::~H5MDTrajectoryReader() {
  delete particle_group_;
  delete atom_position_group_;
  delete atom_force_group_;
  delete edges_group_;
  delete ds_atom_position_;
  delete ds_atom_force_;
  delete ds_atom_velocity_;
  delete ds_time_;
  delete ds_step_;
  delete ds_edges_;
  delete[] time_set_;
  delete[] step_set_;
  delete[] species_;
  delete h5file_;
}

bool H5MDTrajectoryReader::Open(const string &file) {
  // Turn on exceptions instead of various messages.
  H5::Exception::dontPrint();

  // Checks if we deal with hdf5 file.
  if (! H5::H5File::isHdf5(file) )
    throw std::ios_base::failure("Error on open trajectory file: " + file);

  h5file_ = new H5::H5File(file, H5F_ACC_RDONLY);

  // Check the version of the file.
  H5::Group *g_h5md = new H5::Group(h5file_->openGroup("h5md"));
  H5::Attribute *at_version = new H5::Attribute(g_h5md->openAttribute("version"));
  int version[2];
  at_version->read(at_version->getDataType(), &version);
  if (version[0] != 1 || (version[0] == 1 && version[1] > 0)) {
    h5file_->close();
    std::cout << "Major version " << version[0] << std::endl;
    std::cout << "Minor version " << version[1] << std::endl;
    throw std::ios_base::failure("Wrong version of H5MD file.");
  }

  // Clean up.
  g_h5md->close();
  delete at_version;
  delete g_h5md;

  // Checks if the particles group exists and what is the number of members.
  particle_group_ = new H5::Group(h5file_->openGroup("particles"));
  if (particle_group_->getNumObjs() == 0) {
    h5file_->close();
    throw std::ios_base::failure("The particles group is empty.");
  }

  has_velocity_ = false;
  has_force_ = false;
  has_static_box_ = true;
  first_frame_ = true;

  return true;
}

void H5MDTrajectoryReader::Close() {
  h5file_->close();
}

bool H5MDTrajectoryReader::FirstFrame(Topology &top) {  // NOLINT const reference
  // First frame has access to topology so here we will check if the particle_group exists.
  if (first_frame_) {
    first_frame_ = false;

    std::string *particle_group_name = new std::string(top.getParticleGroup());
    if(*particle_group_name == "")
      throw std::ios_base::failure("Missing particle group in topology. Please set h5md_particle_group tag with name \
          attribute set to the particle group.");

    particle_group_ = new H5::Group(h5file_->openGroup("particles"));
    try {
      atom_position_group_ = new H5::Group(particle_group_->openGroup(
            *particle_group_name + "/position"));
    } catch (H5::GroupIException not_found_error) {
      std::cout << "Error: The path '/particles/" << *particle_group_name << "/position' not found. Please check the name ";
      std::cout << "of particle group defined in topology" << std::endl;
      throw not_found_error;
    }

    idx_frame_ = -1;

    ds_atom_position_ = new H5::DataSet(atom_position_group_->openDataSet("value"));
    ds_step_ = new H5::DataSet(atom_position_group_->openDataSet("step"));
    ds_time_ = new H5::DataSet(atom_position_group_->openDataSet("time"));

    // Gets number of particles and dimensions.
    H5::DataSpace fs_atom_position_ = ds_atom_position_->getSpace();
    hsize_t dims[3];
    rank_ = fs_atom_position_.getSimpleExtentDims(dims);
    N_particles_ = dims[1];
    variables_ = dims[2];

    if(N_particles_ != top.BeadCount()) {
      std::cout << "Warning: The number of beads (" << N_particles_ << ")";
      std::cout << " in the trajectory is different than defined in the topology (" << top.BeadCount() << ")" << std::endl;
      std::cout << "The number of beads from topology will be used!" << std::endl;
      N_particles_ = top.BeadCount();
    }
    chunk_rows_[0] = 1;
    chunk_rows_[1] = N_particles_;
    chunk_rows_[2] = variables_;

    // Reads the time set and step set.
    H5::DataSet *ds_pos_step = new H5::DataSet(atom_position_group_->openDataSet("step"));
    H5::DataSet *ds_pos_time = new H5::DataSet(atom_position_group_->openDataSet("time"));
    H5::DataSpace dsp_step = H5::DataSpace(ds_pos_step->getSpace());
    H5::DataSpace dsp_time = H5::DataSpace(ds_pos_time->getSpace());
    hsize_t step_dims[3];
    dsp_step.getSimpleExtentDims(step_dims);
    time_set_ = new double[step_dims[0]];
    step_set_ = new int[step_dims[0]];
    ds_pos_step->read(step_set_, H5::PredType::NATIVE_INT, dsp_step, dsp_step);
    ds_pos_time->read(time_set_, H5::PredType::NATIVE_DOUBLE, dsp_time, dsp_time);

    delete ds_pos_step;
    delete ds_pos_time;

    // Reads the box information.
    H5::Group g_box = H5::Group(particle_group_->openGroup(
        *particle_group_name + "/box"));

    H5::Attribute at_box_dimension = H5::Attribute(g_box.openAttribute("dimension"));
    int dimension;
    at_box_dimension.read(at_box_dimension.getDataType(), &dimension);
    if (dimension != 3 ) {
      throw std::ios_base::failure("Wrong dimension " + dimension);
    }
    // TODO: check also attribute periodic
    try {
      edges_group_ = new H5::Group(g_box.openGroup("edges"));
      has_static_box_ = false;
      throw std::runtime_error("The box is not static during the simulation. Is it NVT ensemble?");
    } catch (H5::GroupIException not_found_error) {
      ds_edges_ = new H5::DataSet(g_box.openDataSet("edges"));
      has_static_box_ = true;
    }

    // Gets the force group.
    try {
      std::string force_group_name = *particle_group_name + "/force";
      atom_force_group_ = new H5::Group(particle_group_->openGroup(force_group_name));
      ds_atom_force_ = new H5::DataSet(atom_force_group_->openDataSet("value"));
    } catch (H5::GroupIException not_found) {
      has_force_ = false;
    }
    // Gets the velocity group.
    try {
      std::string velocity_group_name = *particle_group_name + "/velocity";
      atom_velocity_group_ = new H5::Group(particle_group_->openGroup(velocity_group_name));
      ds_atom_velocity_ = new H5::DataSet(atom_velocity_group_->openDataSet("value"));
    } catch (H5::GroupIException not_found) {
      has_velocity_ = false;
    }
    delete particle_group_name;
  }
  NextFrame(top);
  return true;
}

double* H5MDTrajectoryReader::ReadData(H5::DataSet *ds, int row){
  hsize_t offset[3] = {row, 0, 0};
  H5::DataSpace dsp = H5::DataSpace(ds->getSpace());
  dsp.selectHyperslab(H5S_SELECT_SET, chunk_rows_, offset);
  H5::DataSpace mspace1(variables_, chunk_rows_);
  double *data_out = new double[N_particles_ * variables_];
  ds->read(data_out, H5::PredType::NATIVE_DOUBLE, mspace1, dsp);
  return data_out;
}

/// Reading the data.
bool H5MDTrajectoryReader::NextFrame(Topology &top) {  // NOLINT const reference
  // Reads the position row.
  H5::Exception::dontPrint();
  idx_frame_++;
  std::cout << '\r' << "Reading frame: " << idx_frame_;
  std::cout.flush();
  hsize_t offset[3] = {idx_frame_, 0, 0};
  double *positions;
  try {
    positions = ReadData(ds_atom_position_, idx_frame_);
  } catch (H5::DataSetIException not_found) {
      return false;
  }

  double *forces = NULL;
  double *velocities = NULL;
  if (has_velocity_) {
    velocities = ReadData(ds_atom_velocity_, idx_frame_);
  }

  if (has_force_) {
    forces = ReadData(ds_atom_force_, idx_frame_);
  }

  //Process atoms.
  for (int at_idx = 0; at_idx < N_particles_; at_idx++) {
    double x, y, z;
    int array_index = at_idx*variables_;
    x = positions[array_index];
    y = positions[array_index + 1];
    z = positions[array_index + 2];
    // Topology has to be defined in the xml file or in other
    // topology files. The h5md only stores the trajectory data.
    Bead *b = top.getBead(at_idx);
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

  return true;
}

}  // end csg
}  // end votca

