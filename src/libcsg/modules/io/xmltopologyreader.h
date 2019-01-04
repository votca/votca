/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _XMLTOPOLOGYREADER_H
#define	_XMLTOPOLOGYREADER_H

#include <string>
#include <votca/csg/topologyreader.h>
#include <stack>
#include <votca/tools/parsexml.h>
#include <boost/unordered_map.hpp>


namespace votca { namespace csg {

namespace TOOLS = votca::tools;

class BondBead {
 public:
  BondBead(string &line) {
    TOOLS::Tokenizer tok(line, ":");
    std::vector<std::string> tmp_vec;
    tok.ToVector(tmp_vec);
    if (tmp_vec.size() != 2)
      throw runtime_error("Wrong number of elements in bead: " + line);
    molname = tmp_vec[0];
    atname = tmp_vec[1];
    molname.erase(molname.find_last_not_of(" \n\r\t") + 1);
    atname.erase(atname.find_last_not_of(" \n\r\t") + 1);
  }

  std::string molname;
  std::string atname;
};

class XMLBead {
 public:
  XMLBead(string _name, string _type, double _mass=1.0, double _q=0.0):
    name(_name), type(_type), mass(_mass), q(_q) {};
  XMLBead() {};

  int pid;
  string name;
  string type;
  double mass;
  double q;

};

class XMLMolecule {
 public:
  XMLMolecule(string _name, int _nmols): name(_name), nmols(_nmols) {}
  string name;
  int nmols;
  int pid;
  std::vector<XMLBead*> beads;
  std::map<std::string, XMLBead*> name2beads;
  Molecule *mi;
};

/**
 *  Reads in an xml topology
 *
 * \todo this is a sloppy implementation using expat, is just reads attributes
 * \todo should be extended to also read beads, ...
 * 
*/
class XMLTopologyReader
   : public TopologyReader
{
public:
    /// read a topology file
    bool ReadTopology(std::string file, Topology &top);
    ~XMLTopologyReader();
private:
    typedef boost::unordered_multimap<std::string, XMLMolecule*> MoleculesMap;

    void ReadTopolFile(std::string file);

    void ParseRoot(Property &el);
    void ParseMolecules(Property &el);
    void ParseBeadTypes(Property &el);
    void ParseBonded(Property &el);
    void ParseBox(Property &p);
    void ParseMolecule(Property &p, std::string molname, int nbeads, int nmols);
    void ParseBond(Property &p);
    void ParseAngle(Property &p);
    void ParseDihedral(Property &p);

private:
    ParseXML _parser;

    Topology *_top;
    MoleculesMap _molecules;
    int _mol_index;
    int _bead_index;

    bool _has_base_topology;
};

}}

#endif	/* _PDBTOPOLOGYREADER_H */

