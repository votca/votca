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

#include <numeric>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/map.h>
#include <votca/csg/topology.h>
#include <votca/tools/eigen.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {
class Molecule;
}  // namespace csg
namespace tools {
class Property;
}  // namespace tools
}  // namespace votca

namespace votca {
namespace csg {
using namespace tools;
using namespace std;

Map::~Map() {

  for (auto &_map : _maps) {
    delete _map;
  }
  _maps.clear();
}

void Map::Apply(const BoundaryCondition & bc) {
  for (auto &_map : _maps) {
    _map->Apply(bc);
  }
}

void Map_Sphere::Initialize(Molecule *in, Bead *out, Property *opts_bead,
                            Property *opts_map) {
  BeadMap::Initialize(in, out, opts_bead, opts_map);

  vector<string> beads;
  vector<double> weights;
  vector<double> fweights;

  // get the beads
  string s(_opts_bead->get("beads").value());
  Tokenizer tok_beads(s, " \n\t");
  tok_beads.ToVector(beads);

  // get vector of weights
  Tokenizer tok_weights(_opts_map->get("weights").value(), " \n\t");
  tok_weights.ConvertToVector<double>(weights);

  // check weather weights and # beads matches
  if (beads.size() != weights.size()) {
    throw runtime_error(
        string("number of subbeads in " + opts_bead->get("name").as<string>() +
               " and number of weights in map " +
               opts_map->get("name").as<string>() + " do not match"));
  }

  // normalize the weights
  double norm = 1. / std::accumulate(weights.begin(), weights.end(), 0.);

  transform(weights.begin(), weights.end(), weights.begin(),
            bind2nd(multiplies<double>(), norm));
  // get the d vector if exists or initialize same as weights
  vector<double> d;
  if (_opts_map->exists("d")) {
    Tokenizer tok_weights2(_opts_map->get("d").value(), " \n\t");
    tok_weights2.ConvertToVector(d);
    // normalize d coefficients
    norm = 1. / std::accumulate(d.begin(), d.end(), 0.);
    transform(d.begin(), d.end(), d.begin(),
              bind2nd(multiplies<double>(), norm));
  } else {
    // initialize force-weights with weights
    d.resize(weights.size());
    copy(weights.begin(), weights.end(), d.begin());
  }

  // check weather number of d coeffs is correct
  if (beads.size() != d.size()) {
    throw runtime_error(
        string("number of subbeads in " + opts_bead->get("name").as<string>() +
               " and number of d-coefficients in map " +
               opts_map->get("name").as<string>() + " do not match"));
  }

  fweights.resize(weights.size());
  // calculate force weights by d_i/w_i
  for (size_t i = 0; i < weights.size(); ++i) {
    if (weights[i] == 0 && d[i] != 0) {
      throw runtime_error(
          "A d coefficient is nonzero while weights is zero in mapping " +
          opts_map->get("name").as<string>());
    }
    if (weights[i] != 0) {
      fweights[i] = d[i] / weights[i];
    } else {
      fweights[i] = 0;
    }
  }

  for (size_t i = 0; i < beads.size(); ++i) {
    Index iin = in->getBeadByName(beads[i]);
    if (iin < 0) {
      throw std::runtime_error(
          string("mapping error: molecule " + beads[i] + " does not exist"));
    }
    AddElem(in->getBead(iin), weights[i], fweights[i]);
  }
}

void Map_Sphere::Apply(const BoundaryCondition & bc) {

  bool bPos, bVel, bF;
  bPos = bVel = bF = false;
  _out->ClearParentBeads();

  // the following is needed for pbc treatment
  double max_dist = 0.5 * bc.getShortestBoxDimension();
  Eigen::Vector3d r0 = Eigen::Vector3d::Zero();
  string name0;
  Index id0 = 0;
  if (_matrix.size() > 0) {
    if (_matrix.front()._in->HasPos()) {
      r0 = _matrix.front()._in->getPos();
      name0 = _matrix.front()._in->getName();
      id0 = _matrix.front()._in->getId();
    }
  }

  double M = 0;
  Eigen::Vector3d cg = Eigen::Vector3d::Zero();
  Eigen::Vector3d f = Eigen::Vector3d::Zero();
  Eigen::Vector3d vel = Eigen::Vector3d::Zero();

  for (auto &iter : _matrix) {
    Bead *bead = iter._in;
    _out->AddParentBead(bead->getId());
    M += bead->getMass();
    if (bead->HasPos()) {
      Eigen::Vector3d r = bc.BCShortestConnection(r0, bead->getPos());
      if (r.norm() > max_dist) {
        cout << r0 << " " << bead->getPos() << endl;
        throw std::runtime_error(
            "coarse-grained bead is bigger than half the box \n (atoms " +
            name0 + " (id " + boost::lexical_cast<string>(id0 + 1) + ")" +
            ", " + bead->getName() + " (id " +
            boost::lexical_cast<string>(bead->getId() + 1) + ")" +
            +" , molecule " +
            boost::lexical_cast<string>(bead->getMoleculeId() + 1) + ")");
      }
      cg += iter._weight * (r + r0);
      bPos = true;
    }
    if (bead->HasVel()) {
      vel += iter._weight * bead->getVel();
      bVel = true;
    }
    if (bead->HasF()) {
      f += iter._force_weight * bead->getF();
      bF = true;
    }
  }
  _out->setMass(M);
  if (bPos) {
    _out->setPos(cg);
  }
  if (bVel) {
    _out->setVel(vel);
  }
  if (bF) {
    _out->setF(f);
  }
}

/// \todo implement this function
void Map_Ellipsoid::Apply(const BoundaryConditions & bc) {

  bool bPos, bVel, bF;
  bPos = bVel = bF = false;

  // the following is needed for pbc treatment
  double max_dist = 0.5 * bc.getShortestBoxDimension();
  Eigen::Vector3d r0 = Eigen::Vector3d::Zero();
  if (_matrix.size() > 0) {
    if (_matrix.front()._in->HasPos()) {
      r0 = _matrix.front()._in->getPos();
    }
  }
  Eigen::Vector3d cg = Eigen::Vector3d::Zero();
  Eigen::Vector3d c = Eigen::Vector3d::Zero();
  Eigen::Vector3d f = Eigen::Vector3d::Zero();
  Eigen::Vector3d vel = Eigen::Vector3d::Zero();

  Index n;
  n = 0;
  _out->ClearParentBeads();
  for (auto &iter : _matrix) {
    Bead *bead = iter._in;
    _out->AddParentBead(bead->getId());
    if (bead->HasPos()) {
      Eigen::Vector3d r = bc.BCShortestConnection(r0, bead->getPos());
      if (r.norm() > max_dist) {
        throw std::runtime_error(
            "coarse-grained bead is bigger than half the box");
      }
      cg += iter._weight * (r + r0);
      bPos = true;
    }
    if (bead->HasVel() == true) {
      vel += iter._weight * bead->getVel();
      bVel = true;
    }
    if (bead->HasF()) {
      /// \todo fix me, right calculation should be F_i = m_cg / sum(w_i) *
      /// sum(w_i/m_i*F_i)
      // f += (*iter)._weight * _in->getBeadF((*iter)._in);
      f += iter._force_weight * bead->getF();
      bF = true;
    }

    if (iter._weight > 0 && bead->HasPos()) {
      c += bead->getPos();
      n++;
    }
  }

  if (bPos) {
    _out->setPos(cg);
  }
  if (bVel) {
    _out->setVel(vel);
  }
  if (bF) {
    _out->setF(f);
  }

  if (!_matrix[0]._in->HasPos()) {
    _out->setU(Eigen::Vector3d::UnitX());
    _out->setV(Eigen::Vector3d::UnitY());
    _out->setW(Eigen::Vector3d::UnitZ());
    return;
  }

  Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
  // calculate the tensor of gyration
  c = c / (double)n;
  for (auto &iter : _matrix) {
    if (iter._weight == 0) {
      continue;
    }
    Bead *bead = iter._in;
    Eigen::Vector3d v = bead->getPos() - c;
    // v = vec(1, 0.5, 0) * 0.*(drand48()-0.5)
    //    + vec(0.5, -1, 0) * (drand48()-0.5)
    //    + vec(0, 0, 1) * (drand48()-0.5);

    // Normalize the tensor with 1/number_of_atoms_per_bead
    m += v * v.transpose() / (double)_matrix.size();
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig;
  eig.computeDirect(m);

  Eigen::Vector3d u = eig.eigenvectors().col(0);
  Eigen::Vector3d v = _matrix[1]._in->getPos() - _matrix[0]._in->getPos();
  v.normalize();

  _out->setV(v);

  Eigen::Vector3d w = _matrix[2]._in->getPos() - _matrix[0]._in->getPos();
  w.normalize();

  if (v.cross(w).dot(u) < 0) {
    u = -u;
  }
  _out->setU(u);

  // write out w
  w = u.cross(v);
  w.normalize();
  _out->setW(w);
}

}  // namespace csg
}  // namespace votca
