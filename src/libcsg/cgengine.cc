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

#include <fstream>
#include <votca/csg/cgengine.h>
#include <votca/csg/version.h>
#include <votca/tools/tokenizer.h>

namespace votca {
namespace csg {

using namespace std;

namespace po = boost::program_options;

CGEngine::CGEngine() {}

CGEngine::~CGEngine() {
  map<string, CGMoleculeDef *>::iterator i;
  for (i = _molecule_defs.begin(); i != _molecule_defs.end(); ++i)
    delete (*i).second;
  _molecule_defs.clear();
}

/**
    \todo melts with different molecules
*/
TopologyMap *CGEngine::CreateCGTopology(Topology &in, Topology &out) {
  MoleculeContainer &         mols = in.Molecules();
  MoleculeContainer::iterator iter;
  TopologyMap *               m = new TopologyMap(&in, &out);
  for (iter = mols.begin(); iter != mols.end(); ++iter) {
    Molecule *mol = *iter;
    if (IsIgnored(mol->getName())) continue;
    CGMoleculeDef *def = getMoleculeDef(mol->getName());
    if (!def) {
      cout << "--------------------------------------\n"
           << "WARNING: unknown molecule \"" << mol->getName() << "\" with id "
           << mol->getId() << " in topology" << endl
           << "molecule will not be mapped to CG representation\n"
           << "Check weather a mapping file for all molecule exists, was "
              "specified in --cg "
           << "separated by ; and the ident tag in xml-file matches the "
              "molecule name\n"
           << "--------------------------------------\n";
      continue;
    }
    Molecule *mcg = def->CreateMolecule(out);
    Map *     map = def->CreateMap(*mol, *mcg);
    m->AddMoleculeMap(map);
  }
  out.RebuildExclusions();
  return m;
}

void CGEngine::LoadMoleculeType(string filename) {
  Tokenizer           tok(filename, ";");
  Tokenizer::iterator iter;

  for (iter = tok.begin(); iter != tok.end(); ++iter) {
    CGMoleculeDef *def  = new CGMoleculeDef();
    string         file = *iter;
    boost::trim(file);
    def->Load(file);
    _molecule_defs[def->getIdent()] = def;
  }
}

}  // namespace csg
}  // namespace votca
