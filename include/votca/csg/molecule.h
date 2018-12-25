/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_MOLECULE_H
#define	_VOTCA_CSG_MOLECULE_H

#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include "topologyitem.h"
#include "bead.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

class Interaction;


/**
    \brief information about molecules

    The Molecule class stores which beads belong to a molecule.
    The organization of beads into molecules is needed for the CG mapping.

    \todo sort atoms in molecule

*/
class Molecule : public TopologyItem
{
public:            
    /// get the molecule ID
    int getId() const { return _id; }
    
    /// get the name of the molecule
    const string &getName() const { return _name; }
    
    /// set the name of the molecule
    void setName(const string &name) {  _name=name; }
    
    /// Add a bead to the molecule
    void AddBead(Bead *bead, const string &name);
    /// get the id of a bead in the molecule
    Bead *getBead(int bead) { return _beads[bead]; }
    int getBeadId(int bead) { return _beads[bead]->getId(); }
    int getBeadIdByName(const string &name);
    
    /// get the number of beads in the molecule
    int BeadCount() const { return _beads.size(); }
    
    /// find a bead by it's name
    int getBeadByName(const string &name);
    string getBeadName(int bead) {return _bead_names[bead]; }

    /// Add an interaction to the molecule
    void AddInteraction(Interaction *ic) { _interactions.push_back(ic);
        }

    vector<Interaction *> Interactions() { return _interactions; }

    template<typename T>
    void setUserData(T *userdata) { _userdata = (void*)userdata; }

    template<typename T>
    T *getUserData() { return (T *)_userdata; }
    
private:
    // maps a name to a bead id
    map<string, int> _beadmap;
   vector<Interaction*> _interactions;
     
    // id of the molecules
    int _id;
    
    // name of the molecule
    string _name;
    // the beads in the molecule
    vector<Bead *> _beads;
    vector<string> _bead_names;

    void *_userdata;
    
    /// constructor
    Molecule(Topology *parent, int id, string name)
        : TopologyItem(parent), _id(id), _name(name)
    {}

    friend class Topology;
};

inline int Molecule::getBeadIdByName(const string &name)
{
    int i = getBeadByName(name);
    if(i<0)
        return i;
    return _beads[i]->getId();
}

}}

#endif	/* _VOTCA_CSG_MOLECULE_H */

