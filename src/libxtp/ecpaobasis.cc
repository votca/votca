/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// Standard includes
#include <vector>

// Local VOTCA includes
#include "votca/xtp/checkpointtable.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/ecpbasisset.h"
#include "votca/xtp/qmmolecule.h"
namespace votca {
namespace xtp {

Index ECPAOBasis::getMaxL() const {
  int maxL = 0;
  for (const auto& potential : aopotentials_) {
    maxL = std::max(potential.getL(), maxL);
  }
  return Index(maxL);
}

void ECPAOBasis::AddECPChargeToMolecule(QMMolecule& mol) const {
  for (Index i = 0; i < mol.size(); i++) {
    mol[i].ecpcharge_ = ncore_perAtom_[i];
  }
}

void ECPAOBasis::clear() {
  name_ = "";
  aopotentials_.clear();
  ncore_perAtom_.clear();
}

void ECPAOBasis::UpdatePotentialPositions(const QMMolecule& mol) {
  for (libecpint::ECP& potential : aopotentials_) {
    const Eigen::Vector3d pos = mol[potential.atom_id].getPos();
    potential.center_ = {pos.x(), pos.y(), pos.z()};
  }
}

void ECPAOBasis::add(const ECPAOBasis& other) {
  Index atomindex_offset = Index(ncore_perAtom_.size());
  for (libecpint::ECP potential : other) {
    potential.atom_id += int(atomindex_offset);
    aopotentials_.push_back(potential);
  }

  ncore_perAtom_.insert(ncore_perAtom_.end(), other.ncore_perAtom_.begin(),
                        other.ncore_perAtom_.end());
}

std::vector<std::string> ECPAOBasis::Fill(const ECPBasisSet& bs,
                                          QMMolecule& atoms) {
  aopotentials_.clear();
  ncore_perAtom_.clear();
  name_ = bs.Name();
  std::vector<std::string> non_ecp_elements;
  for (QMAtom& atom : atoms) {
    std::string name = atom.getElement();

    bool element_exists = true;

    try {
      bs.getElement(name);
    } catch (std::runtime_error& error) {
      element_exists = false;

      if (std::find(non_ecp_elements.begin(), non_ecp_elements.end(), name) !=
          non_ecp_elements.end()) {
        non_ecp_elements.push_back(name);
      }
    }

    if (element_exists) {
      const ECPElement& element = bs.getElement(name);
      ncore_perAtom_.push_back(element.getNcore());

      libecpint::ECP potential(atom.getPos().data());
      for (const ECPShell& shell : element) {
        for (const auto& gaussian : shell) {
          potential.addPrimitive(int(gaussian.power_), int(shell.getL()),
                                 gaussian.decay_, gaussian.contraction_);
        }
      }
      aopotentials_.push_back(potential);
      aopotentials_.back().atom_id =
          int(atom.getId());  // add atom id here because copyconstructor  of
                              // libecpint::ECP is broken
    } else {
      ncore_perAtom_.push_back(0);
    }
  }

  AddECPChargeToMolecule(atoms);
  return non_ecp_elements;
}

// This class is purely for hdf5 input output, because the libecpint classes do
// not have the member functions required.
class PotentialIO {

 public:
  struct data {
    Index atomid;
    double x;
    double y;
    double z;
    Index power;
    Index l;
    double coeff;
    double decay;
  };

  static void SetupCptTable(CptTable& table) {
    table.addCol<Index>("atomid", HOFFSET(data, atomid));
    table.addCol<double>("posX", HOFFSET(data, x));
    table.addCol<double>("posY", HOFFSET(data, y));
    table.addCol<double>("posZ", HOFFSET(data, z));
    table.addCol<Index>("power", HOFFSET(data, power));
    table.addCol<Index>("L", HOFFSET(data, l));
    table.addCol<double>("coeff", HOFFSET(data, coeff));
    table.addCol<double>("decay", HOFFSET(data, decay));
  }
};

void ECPAOBasis::WriteToCpt(CheckpointWriter& w) const {
  w(name_, "name");

  w(ncore_perAtom_, "atomic ecp charges");
  Index numofprimitives = 0;
  for (const auto& potential : aopotentials_) {
    numofprimitives += potential.getN();
  }

  CptTable table = w.openTable<PotentialIO>("Potentials", numofprimitives);

  std::vector<PotentialIO::data> dataVec;
  dataVec.reserve(numofprimitives);
  for (const auto& potential : aopotentials_) {
    for (const auto& contrib : potential.gaussians) {
      PotentialIO::data d;
      d.l = Index(contrib.l);
      d.power = Index(contrib.n);
      d.coeff = contrib.d;
      d.decay = contrib.a;
      d.atomid = Index(potential.atom_id);
      d.x = potential.center_[0];
      d.y = potential.center_[1];
      d.z = potential.center_[2];
      dataVec.push_back(d);
    }
  }

  table.write(dataVec);
}

void ECPAOBasis::ReadFromCpt(CheckpointReader& r) {
  clear();
  r(name_, "name");
  r(ncore_perAtom_, "atomic ecp charges");

  CptTable table = r.openTable<PotentialIO>("Potentials");
  std::vector<PotentialIO::data> dataVec(table.numRows());
  table.read(dataVec);
  Index atomindex = -1;
  for (const PotentialIO::data& d : dataVec) {
    if (d.atomid > atomindex) {
      Eigen::Vector3d pos(d.x, d.y, d.z);
      aopotentials_.push_back(libecpint::ECP(pos.data()));
      aopotentials_.back().atom_id = int(d.atomid);
      atomindex = d.atomid;
    }
    // +2 because constructor of libecpint::primitve always subtracts 2
    aopotentials_.back().addPrimitive(int(d.power) + 2, int(d.l), d.decay,
                                      d.coeff);
  }
}

std::ostream& operator<<(std::ostream& out, const libecpint::ECP& potential) {
  out << " AtomId: " << potential.atom_id << " Components: " << potential.getN()
      << "\n";
  for (const auto& gaussian : potential.gaussians) {
    out << " L: " << gaussian.l;
    out << " Gaussian Decay: " << gaussian.a;
    out << " Power: " << gaussian.n;
    out << " Contractions: " << gaussian.d << "\n";
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const ECPAOBasis& ecp) {

  out << "Name:" << ecp.Name() << "\n";
  out << " Potentials:" << ecp.aopotentials_.size() << "\n";
  for (const auto& potential : ecp) {
    out << potential;
  }
  out << "\n"
      << " Atomcharges:";
  for (Index i = 0; i < Index(ecp.ncore_perAtom_.size()); i++) {
    out << i << ":" << ecp.ncore_perAtom_[i] << " ";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
