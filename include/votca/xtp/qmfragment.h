/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include <boost/lexical_cast.hpp>
#include <limits>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/bse_population.h>

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
  QMFragment(std::string name, int id, std::string atoms)
      : _name(name), _id(id) {
    FillAtomIndices(atoms);
  }

  QMFragment(){};

  QMFragment(CheckpointReader& r) { ReadFromCpt(r); }

  void setName(std::string name) { _name = name; }
  void setId(int id) { _id = id; }
  void getName() const { return _name; }
  void getId() const { return _id; }
  void FillFromString(std::string atoms) { FillAtomIndices(atoms); }

  const T& value() const { return _value; }

  T& value() { return _value; }
  const std::string& name() const { return _name; }

  int size() const { return _atomindices.size(); }

  const std::vector<int>& getIndices() const { return _atomindices; }

  double ExtractFromVector(const Eigen::VectorXd& atomentries) const {
    double result = 0;
    for (int index : _atomindices) {
      result += atomentries(index);
    }
    return result;
  }

  typename std::vector<int>::const_iterator begin() const {
    return _atomindices.begin();
  }
  typename std::vector<int>::const_iterator end() const {
    return _atomindices.end();
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const QMFragment& fragment) {
    out << "Fragment name:" << fragment._name << " id:" << fragment._id << "\n";
    out << "AtomIndices[" << fragment.size() << "]:";
    for (int id : fragment._atomindices) {
      out << id << " ";
    }
    out << "\nValue:" << fragment._value;
    out << "\n";
    return out;
  };

  void WriteToCpt(CheckpointWriter& w) const {
    w(_atomindices, "indices");
    w(_name, "name");
    w(_id, "id");
    WriteValue(w);
  }

  void ReadFromCpt(CheckpointReader& r) {
    r(_atomindices, "indices");
    r(_name, "name");
    r(_id, "id");
    ReadValue(r);
  }

 private:
  void WriteValue(CheckpointWriter& w) const;
  void ReadValue(CheckpointReader& r);

  void FillAtomIndices(const std::string& atoms) {
    tools::Tokenizer tok(atoms, " ,\n\t");
    std::vector<std::string> results;
    tok.ToVector(results);
    const std::string delimiter = "...";
    for (std::string s : results) {
      if (s.find(delimiter) != std::string::npos) {
        int start = std::stoi(s.substr(0, s.find(delimiter)));
        int stop =
            std::stoi(s.erase(0, s.find(delimiter) + delimiter.length()));
        for (int i = start; i <= stop; i++) {
          _atomindices.push_back(i);
        }
      } else {
        _atomindices.push_back(std::stoi(s));
      }
    }
  }

  std::vector<int> _atomindices;
  std::string _name = "";
  int _id = -1;
  T _value;
};

template <class T>
inline void QMFragment<T>::ReadValue(CheckpointReader& r) {
  r(_value, "value");
}
template <class T>
inline void QMFragment<T>::WriteValue(CheckpointWriter& w) const {
  w(_value, "value");
}

template <>
inline void QMFragment<BSE_Population>::ReadValue(CheckpointReader& r) {
  CheckpointReader rr = r.openChild("BSE_pop");
  _value.ReadFromCpt(rr);
}

template <>
inline void QMFragment<BSE_Population>::WriteValue(CheckpointWriter& w) const {
  CheckpointWriter ww = w.openChild("BSE_pop");
  _value.WriteToCpt(ww);
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOMCONTAINER_H
