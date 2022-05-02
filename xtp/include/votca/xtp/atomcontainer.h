

/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_ATOMCONTAINER_H
#define VOTCA_XTP_ATOMCONTAINER_H

// Standard includes
#include <limits>
#include <typeinfo>

// VOTCA includes
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "checkpoint.h"
#include "eigen.h"

/**
 * \brief Basic Container for QMAtoms,PolarSites and Atoms
 *
 *
 *
 */

namespace votca {
namespace xtp {

template <class T>
class AtomContainer {
 public:
  AtomContainer(std::string type, Index id) : type_(type), id_(id){};

  using Atom_Type = T;

  AtomContainer(CheckpointReader& r) { this->ReadFromCpt(r); }
  virtual ~AtomContainer() = default;

  using iterator = typename std::vector<T>::iterator;

  const std::string& getType() const { return type_; }

  void setType(std::string type) { type_ = type; }

  void clearAtoms() { atomlist_.clear(); }

  Index getId() const { return id_; }

  Index size() const { return atomlist_.size(); }

  void push_back(const T& atom) {
    atomlist_.push_back(atom);
    calcPos();
  }
  void push_back(T&& atom) {
    atomlist_.push_back(atom);
    calcPos();
  }

  const T& at(Index index) const { return atomlist_.at(index); }
  T& at(Index index) { return atomlist_.at(index); }

  const T& operator[](Index index) const { return atomlist_[index]; }
  T& operator[](Index index) { return atomlist_[index]; }

  typename std::vector<T>::iterator begin() { return atomlist_.begin(); }
  typename std::vector<T>::iterator end() { return atomlist_.end(); }

  typename std::vector<T>::const_iterator begin() const {
    return atomlist_.begin();
  }
  typename std::vector<T>::const_iterator end() const {
    return atomlist_.end();
  }

  const Eigen::Vector3d& getPos() const { return pos_; }

  // calculates the lowest and highest point in the cube, sorrounding the
  // molecule
  std::pair<Eigen::Vector3d, Eigen::Vector3d> CalcSpatialMinMax() const {
    std::pair<Eigen::Vector3d, Eigen::Vector3d> result;
    Eigen::Vector3d min =
        std::numeric_limits<double>::max() * Eigen::Vector3d::Ones();
    Eigen::Vector3d max =
        std::numeric_limits<double>::min() * Eigen::Vector3d::Ones();
    for (const T& atom : atomlist_) {
      const Eigen::Vector3d& pos = atom.getPos();
      if (pos.x() < min.x()) {
        min.x() = pos.x();
      }
      if (pos.x() > max.x()) {
        max.x() = pos.x();
      }
      if (pos.y() < min.y()) {
        min.y() = pos.y();
      }
      if (pos.y() > max.y()) {
        max.y() = pos.y();
      }
      if (pos.z() < min.z()) {
        min.z() = pos.z();
      }
      if (pos.z() > max.z()) {
        max.z() = pos.z();
      }
    }
    result.first = min;
    result.second = max;
    return result;
  }

  std::vector<std::string> FindUniqueElements() const {
    std::vector<std::string> result;
    for (const T& atom : atomlist_) {
      if (std::find(result.begin(), result.end(), atom.getElement()) ==
          result.end()) {
        result.push_back(atom.getElement());
      }
    }
    return result;
  }

  void Translate(const Eigen::Vector3d& shift) {
    for (T& atom : atomlist_) {
      atom.Translate(shift);
    }
    pos_ += shift;
  }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos) {
    for (T& atom : atomlist_) {
      atom.Rotate(R, ref_pos);
    }
    calcPos();
  }

  virtual void WriteToCpt(CheckpointWriter& w) const {
    w(type_, "type");
    w(id_, "id");
    w(int(atomlist_.size()), "size");
    T element(0, "H", Eigen::Vector3d::Zero());
    CptTable table = w.openTable<T>(element.identify() + "s", atomlist_.size());
    std::vector<typename T::data> dataVec(atomlist_.size());
    for (std::size_t i = 0; i < atomlist_.size(); ++i) {
      atomlist_[i].WriteData(dataVec[i]);
    }

    table.write(dataVec);
  }
  virtual void ReadFromCpt(CheckpointReader& r) {
    r(type_, "type");
    r(id_, "id");
    Index size = 0;
    r(size, "size");
    if (size == 0) {
      return;
    }
    T element(0, "H", Eigen::Vector3d::Zero());  // dummy element to get
                                                 // .identify for type
    CptTable table = r.openTable<T>(element.identify() + "s");
    atomlist_.clear();
    atomlist_.reserve(table.numRows());
    std::vector<typename T::data> dataVec(table.numRows());
    table.read(dataVec);
    for (std::size_t i = 0; i < table.numRows(); ++i) {
      atomlist_.push_back(T(dataVec[i]));
    }
    calcPos();
  }

  void calcPos() {
    tools::Elements element;
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    double totalmass = 0.0;
    for (const T& atom : atomlist_) {
      double mass = element.getMass(atom.getElement());
      totalmass += mass;
      pos += mass * atom.getPos();
    }
    pos_ = pos / totalmass;
  }

 protected:
  std::vector<T> atomlist_;
  std::string type_;
  Index id_;

 private:
  Eigen::Vector3d pos_ = Eigen::Vector3d::Zero();
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOMCONTAINER_H
