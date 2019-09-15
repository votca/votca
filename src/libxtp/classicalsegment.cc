/*
 *            Copyright 2016 The MUSCET Development Team
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

#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/atomcontainer.h"
#include <boost/format.hpp>
#include <votca/tools/elements.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace xtp {

// MPS files have a weird format positions can be in bohr or angstroem,
// multipoles are in q*bohr^k, with k rank of multipole and polarisabilities are
// in angstroem^3
template <class T>
void ClassicalSegment<T>::LoadFromFile(std::string filename) {

  std::string line;
  std::ifstream intt;
  intt.open(filename);
  double unit_conversion = tools::conv::ang2bohr;

  int readinmultipoles = 0;
  int numberofmultipoles = 0;
  Vector9d multipoles = Vector9d::Zero();
  int rank = 0;

  if (!intt.is_open()) {
    throw std::runtime_error("File:" + filename + " could not be opened");
  }
  while (intt.good()) {

    std::getline(intt, line);
    std::vector<std::string> split;
    tools::Tokenizer toker(line, " \t");
    toker.ToVector(split);

    if (!split.size() || split[0] == "!" || split[0].substr(0, 1) == "!") {
      continue;
    }

    // ! Interesting information here, e.g.
    // ! DCV2T opt
    // ! SP        RB3LYP          6-311+G(d,p)
    // Units bohr
    //
    // C          -4.2414603400   -3.8124751600    0.0017575736    Rank  2
    //  -0.3853409355
    //  -0.0002321905   0.2401559510   0.6602334308
    //  -0.7220625314   0.0004894995  -0.0003833545   0.4526409813  -0.50937399
    //  P 1.75

    // Units used
    if (split[0] == "Units") {
      std::string units = split[1];
      if (units != "bohr" && units != "angstrom") {
        throw std::runtime_error("Unit " + units + " in file " + filename +
                                 " not supported.");
      }
      if (units == "bohr") {
        unit_conversion = 1.0;
      }
    } else if (split.size() == 6) {
      // element,  position,  rank limit convert to bohr
      std::string name = split[0];
      Eigen::Vector3d pos;
      int id = this->_atomlist.size();
      pos[0] = boost::lexical_cast<double>(split[1]);
      pos[1] = boost::lexical_cast<double>(split[2]);
      pos[2] = boost::lexical_cast<double>(split[3]);
      rank = boost::lexical_cast<int>(split[5]);
      numberofmultipoles = (rank + 1) * (rank + 1);
      multipoles = Vector9d::Zero();
      pos *= unit_conversion;
      this->_atomlist.push_back(T(id, name, pos));
    }
    // 'P', dipole polarizability
    else if (split[0] == "P") {
      Eigen::Matrix3d p1;
      if (split.size() == 7) {
        double pxx = boost::lexical_cast<double>(split[1]);
        double pxy = boost::lexical_cast<double>(split[2]);
        double pxz = boost::lexical_cast<double>(split[3]);
        double pyy = boost::lexical_cast<double>(split[4]);
        double pyz = boost::lexical_cast<double>(split[5]);
        double pzz = boost::lexical_cast<double>(split[6]);
        p1 << pxx, pxy, pxz, pxy, pyy, pyz, pxz, pyz, pzz;
      } else if (split.size() == 2) {
        double pxx = boost::lexical_cast<double>(split[1]);
        p1 = pxx * Eigen::Matrix3d::Identity();
      } else {
        throw std::runtime_error("Invalid line in " + filename + ": " + line);
      }
      double unit_conversion_3 = std::pow(tools::conv::ang2bohr, 3);
      p1 = p1 * unit_conversion_3;
      this->_atomlist.back().setPolarisation(p1);
    }
    // Multipole lines
    else {
      // stay in bohr
      for (unsigned i = 0; i < split.size(); i++) {
        double qXYZ = boost::lexical_cast<double>(split[i]);
        if (multipoles.size() < readinmultipoles) {
          throw std::runtime_error("ReadMpsFile: File" + filename +
                                   "is not properly formatted");
        }
        multipoles(readinmultipoles) = qXYZ;
        readinmultipoles++;
      }
      if (readinmultipoles == numberofmultipoles) {
        Eigen::Vector3d temp_dipoles = multipoles.segment<3>(1);
        // mps format for dipoles is z x y
        // we need x y z
        multipoles(1) = temp_dipoles(1);
        multipoles(2) = temp_dipoles(2);
        multipoles(3) = temp_dipoles(0);
        this->_atomlist.back().setMultipole(multipoles, rank);
        readinmultipoles = 0;
      }
    }
  }
  this->calcPos();
}

template <class T>
double ClassicalSegment<T>::CalcTotalQ() const {
  double Q = 0;
  for (const T& site : this->_atomlist) {
    Q += site.getCharge();
  }
  return Q;
}

template <class T>
Eigen::Vector3d ClassicalSegment<T>::CalcDipole() const {
  Eigen::Vector3d dipole = Eigen::Vector3d::Zero();

  Eigen::Vector3d CoM = this->getPos();
  for (const T& site : this->_atomlist) {
    dipole += (site.getPos() - CoM) * site.getCharge();
    dipole += site.getDipole();
  }
  return dipole;
}

template <class T>
void ClassicalSegment<T>::WriteMPS(std::string filename,
                                   std::string header) const {

  std::ofstream ofs;
  ofs.open(filename, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + filename);
  }

  ofs << (boost::format("! GENERATED BY VOTCA::XTP::%1$s\n") % header);
  ofs << (boost::format("! N=%2$d Q[e]=%1$+1.7f\n") % CalcTotalQ() %
          this->size());
  ofs << boost::format("Units angstrom\n");

  for (const T& site : this->_atomlist) {
    ofs << site.WriteMpsLine("angstrom");
  }
  ofs.close();
}

template <class T>
std::string ClassicalSegment<T>::identify() const {
  return "";
}

template <>
std::string ClassicalSegment<PolarSite>::identify() const {
  return "PolarSegment";
}
template <>
std::string ClassicalSegment<StaticSite>::identify() const {
  return "StaticSegment";
}

template class ClassicalSegment<PolarSite>;
template class ClassicalSegment<StaticSite>;

}  // namespace xtp
}  // namespace votca
