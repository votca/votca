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
#ifndef VOTCA_XTP_QMFRAGMENT_H
#define VOTCA_XTP_QMFRAGMENT_H

// Local VOTCA includes
#include "IndexParser.h"
#include "bse_population.h"

/**
 * \brief Container to define fragments of QMmolecules, containing atomindices,
 * no pointers to atoms, it also handles the parsing of strings etc.. Values
 * should have own destructor
 *
 *
 *
 */

namespace votca {
namespace xtp {

template <class T>
class QMFragment {
 public:
  QMFragment(Index id, std::string atoms) : id_(id) { FillAtomIndices(atoms); }

  QMFragment() = default;

  QMFragment(CheckpointReader& r) { ReadFromCpt(r); }

  template <class T2>
  void copy_withoutvalue(const QMFragment<T2>& frag) {
    id_ = frag.getId();
    atomindices_ = frag.getIndices();
  }

  void setId(Index id) { id_ = id; }
  Index getId() const { return id_; }
  void FillFromString(std::string atoms) { FillAtomIndices(atoms); }

  const T& value() const { return value_; }

  T& value() { return value_; }

  Index size() const { return Index(atomindices_.size()); }

  const std::vector<Index>& getIndices() const { return atomindices_; }

  double ExtractFromVector(const Eigen::VectorXd& atomentries) const {
    double result = 0;
    for (Index index : atomindices_) {
      result += atomentries(index);
    }
    return result;
  }

  typename std::vector<Index>::const_iterator begin() const {
    return atomindices_.begin();
  }
  typename std::vector<Index>::const_iterator end() const {
    return atomindices_.end();
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const QMFragment& fragment) {
    out << "Fragment id:" << fragment.id_ << "\n";
    out << "AtomIndices[" << fragment.size() << "]:";
    IndexParser p;
    out << p.CreateIndexString(fragment.atomindices_) << "\n";
    out << "\nValue:" << fragment.value_;
    out << "\n";
    return out;
  };

  void WriteToCpt(CheckpointWriter& w) const {
    w(atomindices_, "indices");
    w(id_, "id");
    WriteValue(w);
  }

  void ReadFromCpt(CheckpointReader& r) {
    r(atomindices_, "indices");
    r(id_, "id");
    ReadValue(r);
  }

 private:
  void WriteValue(CheckpointWriter& w) const;
  void ReadValue(CheckpointReader& r);

  void FillAtomIndices(const std::string& atoms) {
    IndexParser p;
    atomindices_ = p.CreateIndexVector(atoms);
  }

  std::vector<Index> atomindices_;
  Index id_ = -1;
  T value_{};
};

template <class T>
inline void QMFragment<T>::ReadValue(CheckpointReader& r) {
  r(value_, "value");
}
template <class T>
inline void QMFragment<T>::WriteValue(CheckpointWriter& w) const {
  w(value_, "value");
}

template <>
inline void QMFragment<BSE_Population>::ReadValue(CheckpointReader& r) {
  CheckpointReader rr = r.openChild("BSE_pop");
  value_.ReadFromCpt(rr);
}

template <>
inline void QMFragment<BSE_Population>::WriteValue(CheckpointWriter& w) const {
  CheckpointWriter ww = w.openChild("BSE_pop");
  value_.WriteToCpt(ww);
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMFRAGMENT_H
