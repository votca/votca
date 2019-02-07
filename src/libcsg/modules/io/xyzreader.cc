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

#include "xyzreader.h"
#include <boost/lexical_cast.hpp>
#include <vector>
#include <votca/tools/getline.h>

namespace votca {
namespace csg {
    using namespace boost;
using namespace std;

bool XYZReader::ReadTopology(string file, Topology &top) {
   top.Cleanup();

   _fl.open(file.c_str());
  if (!_fl.is_open()) {
        throw std::ios_base::failure("Error on open topology file: " + file);
  }

  ReadFrame<true>(top);

    _fl.close();

    return true;
}

bool XYZReader::Open(const string &file) {
    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open trajectory file: " + file);
    _line = 0;
    return true;
}

void XYZReader::Close() { _fl.close(); }

bool XYZReader::FirstFrame(Topology &top) { return NextFrame(top); }

bool XYZReader::NextFrame(Topology &top) {
  bool success = ReadFrame<false>(top);
  return success;
}

template <bool topology>
bool XYZReader::ReadFrame(Topology &top) {
    string line;
  getline(_fl, line);
  ++_line;
    if(!_fl.eof()) {
        // read the number of atoms
    int natoms = boost::lexical_cast<int>(line);
    if (!topology && natoms != top.BeadCount())
      throw std::runtime_error(
          "number of beads in topology and trajectory differ");

        // the title line
    getline(_fl, line);
    ++_line;

        // read atoms
    for (int i = 0; i < natoms; ++i) {
      getline(_fl, line);
      ++_line;
      if (_fl.eof()) {
                throw std::runtime_error("unexpected end of file in xyz file");
      }

            vector<string> fields;
            Tokenizer tok(line, " ");
            tok.ToVector(fields);

      if (fields.size() != 4) {
                throw std::runtime_error("invalide line " + 
                        boost::lexical_cast<string>(_line) +
                        " in xyz file\n" + line);
      }

            Bead *b;
      if (topology) {
        string bead_type = fields[0];
        if (!top.BeadTypeExist(bead_type)) {
          top.RegisterBeadType(bead_type);
        }
                b = top.CreateBead(1, fields[0]+boost::lexical_cast<string>(i),
                           bead_type, 0, 0, 0);
      } else {
                b = top.getBead(i);
      }

            // convert to nm from A
            b->setPos(vec(
                    stod(fields[1])/10.0,
                    stod(fields[2])/10.0,
                    stod(fields[3])/10.0
                ));

    }  // for(int i=0; i<natoms; ++i)
  }    // if(!_fl.eof())
  return !_fl.eof();
        }

}  // namespace csg
}  // namespace votca
