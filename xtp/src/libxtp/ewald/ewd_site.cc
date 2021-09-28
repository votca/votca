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

// Local VOTCA includes

#include "votca/xtp/segmentmapper.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "votca/xtp/background.h"
#include "votca/xtp/ewd_site.h"
#include "votca/xtp/unitcell.h"

namespace votca {
namespace xtp {

EwdSite::EwdSite(const PolarSite& pol) {
  id_ = pol.getId();
  rank_ = pol.getRank();
  position_ = pol.getPos();
  charge_ = pol.getCharge();
  dipole_static_ = pol.getStaticDipole();
  dipole_induced_ = pol.getInducedDipole();
  polarization_ = pol.getpolarization();
  element_ = pol.getElement();
  // the 1/3 saves factors in further calculations
  // The quadrupole in the Ewald site should by multiplied by 3 to obtain the
  // true quadrupole
  quadrupole_ = (1.0 / 3.0) * pol.CalculateCartesianMultipole();
}

void EwdSite::WriteData(data& d) {
  d.id = id_;
  d.rank = rank_;
  d.posX = position_[0];
  d.posY = position_[1];
  d.posZ = position_[2];

  d.charge = charge_;
  d.dipX = dipole_static_[0];
  d.dipY = dipole_static_[1];
  d.dipZ = dipole_static_[2];
  d.quadXX = quadrupole_(0, 0);
  d.quadXY = quadrupole_(0, 1);
  d.quadXZ = quadrupole_(0, 2);
  d.quadYY = quadrupole_(1, 1);
  d.quadYZ = quadrupole_(1, 2);
  d.quadZZ = quadrupole_(2, 2);

  d.polXX = polarization_(0, 0);
  d.polXY = polarization_(0, 1);
  d.polXZ = polarization_(0, 2);
  d.polYY = polarization_(1, 1);
  d.polYZ = polarization_(1, 2);
  d.polZZ = polarization_(2, 2);

  d.d_x_ind = dipole_induced_[0];
  d.d_y_ind = dipole_induced_[1];
  d.d_z_ind = dipole_induced_[2];

  d.fieldX_stat = field_static_[0];
  d.fieldX_ind = field_induced_[0];
  d.fieldY_stat = field_static_[1];
  d.fieldY_ind = field_induced_[1];
  d.fieldZ_stat = field_static_[2];
  d.fieldZ_ind = field_induced_[2];

  d.element = const_cast<char*>(element_.c_str());
}

void EwdSite::ReadData(const data& d) {
  id_ = d.id;
  rank_ = d.rank;
  position_[0] = d.posX;
  position_[1] = d.posY;
  position_[2] = d.posZ;

  charge_ = d.charge;
  dipole_static_[0] = d.dipX;
  dipole_static_[1] = d.dipY;
  dipole_static_[2] = d.dipZ;
  quadrupole_(0, 0) = d.quadXX;
  quadrupole_(0, 1) = d.quadXY;
  quadrupole_(0, 2) = d.quadXZ;
  quadrupole_(1, 1) = d.quadYY;
  quadrupole_(1, 2) = d.quadYZ;
  quadrupole_(2, 2) = d.quadZZ;
  quadrupole_(2, 0) = d.quadXY;
  quadrupole_(2, 0) = d.quadXZ;
  quadrupole_(2, 1) = d.quadYZ;

  polarization_(0, 0) = d.polXX;
  polarization_(0, 1) = d.polXY;
  polarization_(0, 2) = d.polXZ;
  polarization_(1, 1) = d.polYY;
  polarization_(1, 2) = d.polYZ;
  polarization_(2, 2) = d.polZZ;
  polarization_(2, 0) = d.polXY;
  polarization_(2, 0) = d.polXZ;
  polarization_(2, 1) = d.polYZ;

  dipole_induced_[0] = d.d_x_ind;
  dipole_induced_[1] = d.d_y_ind;
  dipole_induced_[2] = d.d_z_ind;

  field_static_[0] = d.fieldX_stat;
  field_induced_[0] = d.fieldX_ind;
  field_static_[1] = d.fieldY_stat;
  field_induced_[1] = d.fieldY_ind;
  field_static_[2] = d.fieldZ_stat;
  field_induced_[2] = d.fieldZ_ind;
  element_ = std::string(d.element);
}

void EwdSite::SetupCptTable(CptTable& table) {
  table.addCol<Index>("index", HOFFSET(data, id));
  table.addCol<Index>("rank", HOFFSET(data, rank));
  table.addCol<double>("posX", HOFFSET(data, posX));
  table.addCol<double>("posY", HOFFSET(data, posY));
  table.addCol<double>("posZ", HOFFSET(data, posZ));

  table.addCol<double>("charge", HOFFSET(data, charge));
  table.addCol<double>("dipX", HOFFSET(data, dipX));
  table.addCol<double>("dipY", HOFFSET(data, dipY));
  table.addCol<double>("dipZ", HOFFSET(data, dipZ));
  table.addCol<double>("quadXX", HOFFSET(data, quadXX));
  table.addCol<double>("quadXY", HOFFSET(data, quadXY));
  table.addCol<double>("quadXZ", HOFFSET(data, quadXZ));
  table.addCol<double>("quadYY", HOFFSET(data, quadYY));
  table.addCol<double>("quadYZ", HOFFSET(data, quadYZ));
  table.addCol<double>("quadZZ", HOFFSET(data, quadZZ));

  table.addCol<double>("polXX", HOFFSET(data, polXX));
  table.addCol<double>("polXY", HOFFSET(data, polXY));
  table.addCol<double>("polXZ", HOFFSET(data, polXZ));
  table.addCol<double>("polYY", HOFFSET(data, polYY));
  table.addCol<double>("polYZ", HOFFSET(data, polYZ));
  table.addCol<double>("polZZ", HOFFSET(data, polZZ));

  table.addCol<double>("d_x_ind", HOFFSET(data, d_x_ind));
  table.addCol<double>("d_y_ind", HOFFSET(data, d_y_ind));
  table.addCol<double>("d_z_ind", HOFFSET(data, d_z_ind));

  table.addCol<double>("fieldX_stat", HOFFSET(data, fieldX_stat));
  table.addCol<double>("fieldY_stat", HOFFSET(data, fieldY_stat));
  table.addCol<double>("fieldZ_stat", HOFFSET(data, fieldZ_stat));
  table.addCol<double>("fieldX_ind", HOFFSET(data, fieldX_ind));
  table.addCol<double>("fieldY_ind", HOFFSET(data, fieldY_ind));
  table.addCol<double>("fieldZ_ind", HOFFSET(data, fieldZ_ind));
  table.addCol<std::string>("element", HOFFSET(data, element));
}

void EwdSite::induceDirect() {
  dipole_induced_ = -polarization_ * field_static_;
}

}  // namespace xtp
}  // namespace votca