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

#include <boost/format.hpp>
#include <stdio.h>
#include <string>
#include <votca/csg/pdbwriter.h>

namespace votca {
namespace csg {

using namespace std;

void PDBWriter::Open(string file, bool bAppend) {
  if (bAppend) {
    _out.open(file, std::ios_base::app);
  } else {
    _out.open(file);
  }
}

void PDBWriter::WriteHeader(std::string header) {
  if (header.size() < 10 || header.substr(0, 10) != "HEADER    ") {
    _out << "HEADER    ";
  }
  _out << header;
  if (header.back() != '\n') _out << "\n";
}

void PDBWriter::Close() { _out.close(); }

void PDBWriter::Write(Topology *conf) {

  _out << boost::format("MODEL     %1$4d\n") % (conf->getStep() + 1)
       << std::flush;
  ;
  WriteContainer<Topology>(*conf);
  _out << "ENDMDL" << std::endl;
}

void PDBWriter::writeSymmetry(Bead *bead) {
  if (bead->getSymmetry() > 1) {
    Eigen::Vector3d r = 10 * bead->getPos();
    boost::format   beadfrmt(
        "HETATM%1$5d %2$4s %3$3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f\n");
    Eigen::Vector3d ru = 0.1 * bead->getU() + r;

    _out << beadfrmt % (bead->getId() + 1) % 100000  // atom serial number
                % bead->getName()                    // atom name
                % "REU"                              // residue name
                % " "                                // chain identifier 1 char
                % (bead->getResnr() + 1)             // residue sequence number
                % ru.x() % ru.y() % ru.z();          // we skip the charge

    if (bead->getSymmetry() > 2) {
      Eigen::Vector3d rv = 0.1 * bead->getV() + r;
      _out << beadfrmt % (bead->getId() + 1) % 100000  // atom serial number
                  % bead->getName()                    // atom name
                  % "REV"                              // residue name
                  % " "                        // chain identifier 1 char
                  % (bead->getResnr() + 1)     // residue sequence number
                  % rv.x() % rv.y() % rv.z();  // we skip the charge
    }
  }
  return;
}

std::string PDBWriter::getResname(Topology &conf, Bead *bead) {
  if (conf.getResidue(bead->getResnr())) {
    return conf.getResidue(bead->getResnr())->getName();
  } else {
    return "";
  }
}
}  // namespace csg
}  // namespace votca
