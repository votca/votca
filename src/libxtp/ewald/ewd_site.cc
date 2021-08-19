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
#include "backgroundpolarizer.h"
#include "ewd_site.h"
#include "unitcell.h"

namespace votca {
namespace xtp {

EwdSite::EwdSite(const PolarSite& pol) {
  _id = pol.getId();
  _rank = pol.getRank();
  _position = pol.getPos();
  _charge = pol.getCharge();
  _dipole_static = pol.getStaticDipole();
  _dipole_induced = Eigen::Vector3d::Zero();
  _polarization = pol.getpolarization();
  _element = pol.getElement();
  // the 1/3 saves factors in further calculations
  // The quadrupole in the Ewald site should by multiplied by 3 to obtain the
  // true quadrupole
  _quadrupole = (1.0 / 3.0) * pol.CalculateCartesianMultipole();
}

void EwdSite::WriteData(data& d) {
  d.id = _id;
  d.rank = _rank;
  d.posX = _position[0];
  d.posY = _position[1];
  d.posZ = _position[2];

  d.charge = _charge;
  d.dipX = _dipole_static[0];
  d.dipY = _dipole_static[1];
  d.dipZ = _dipole_static[2];
  d.quadXX = _quadrupole(0,0);
  d.quadXY = _quadrupole(0,1);
  d.quadXZ = _quadrupole(0,2);
  d.quadYY = _quadrupole(1,1);
  d.quadYZ = _quadrupole(1,2);
  d.quadZZ = _quadrupole(2,2);

  d.d_x_ind = _dipole_induced[0];
  d.d_y_ind = _dipole_induced[1];
  d.d_z_ind = _dipole_induced[2];

  d.fieldX_stat = _field_static[0];
  d.fieldX_ind = _field_induced[0];
  d.fieldY_stat = _field_static[1];
  d.fieldY_ind = _field_induced[1];
  d.fieldZ_stat = _field_static[2];
  d.fieldZ_ind = _field_induced[2];
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

  table.addCol<double>("d_x_ind", HOFFSET(data, d_x_ind));
  table.addCol<double>("d_y_ind", HOFFSET(data, d_y_ind));
  table.addCol<double>("d_z_ind", HOFFSET(data, d_z_ind));

  table.addCol<double>("fieldX_stat", HOFFSET(data, fieldX_stat));
  table.addCol<double>("fieldY_stat", HOFFSET(data, fieldY_stat));
  table.addCol<double>("fieldZ_stat", HOFFSET(data, fieldZ_stat));
  table.addCol<double>("fieldX_ind", HOFFSET(data, fieldX_ind));
  table.addCol<double>("fieldY_ind", HOFFSET(data, fieldY_ind));
  table.addCol<double>("fieldZ_ind", HOFFSET(data, fieldZ_ind));
}

void EwdSite::induceDirect() { _dipole_induced = -_polarization * _field_static; }

}  // namespace xtp
}  // namespace votca