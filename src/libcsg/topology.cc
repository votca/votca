/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/csg/topology.h>
#include <votca/csg/interaction.h>
#include <votca/tools/rangeparser.h>
#include <stdexcept>

namespace votca { namespace csg {

using namespace std;

Topology::~Topology()
{
    Cleanup();
    if(_bc)
        delete (_bc);
    _bc = NULL;
}

void Topology::Cleanup()
{
    // cleanup beads
    {
        BeadContainer::iterator i;
        for(i=_beads.begin();i<_beads.end();++i)
            delete *i;
        _beads.clear();
    }
    // cleanup molecules
    {
        MoleculeContainer::iterator i;
        for(i=_molecules.begin();i<_molecules.end();++i)
            delete *i;
        _molecules.clear();
    }
    // cleanup residues
    {
        ResidueContainer::iterator i;
        for(i=_residues.begin();i<_residues.end();++i)
            delete (*i);
        _residues.clear();
    }
    // cleanup interactions
    {
        InteractionContainer::iterator i;
        for(i=_interactions.begin();i<_interactions.end();++i)
            delete (*i);
        _interactions.clear();
    }
    // cleanup _bc object
    if(_bc)
        delete (_bc);
    _bc = new OpenBox();
}

/// \todo implement checking, only used in xml topology reader
void Topology::CreateMoleculesByRange(string name, int first, int nbeads, int nmolecules)
{
    Molecule *mol = CreateMolecule(name);
    int beadcount=0;
    int res_offset=0;
   
    BeadContainer::iterator bead;
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
        //xml numbering starts with 1
        if(--first > 0) continue;
	//This is not 100% correct, but let's assume for now that the resnr do increase
	if ( beadcount == 0 ) {
	    res_offset = (*bead)->getResnr();
	}
        stringstream bname;
	bname << (*bead)->getResnr() - res_offset + 1 << ":" << getResidue((*bead)->getResnr())->getName() << ":" << (*bead)->getName();
        mol->AddBead((*bead), bname.str());
        if(++beadcount == nbeads) {
            if(--nmolecules <= 0) break;
            mol = CreateMolecule(name);            
            beadcount = 0;
        }
    }
}

/// \todo clean up CreateMoleculesByResidue!
void Topology::CreateMoleculesByResidue()
{
    // first create a molecule for each residue    
    ResidueContainer::iterator res;    
    for(res=_residues.begin(); res!=_residues.end(); ++res) {
        CreateMolecule((*res)->getName());
    }
    
    // add the beads to the corresponding molecules based on their resid
    BeadContainer::iterator bead;    
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
        //MoleculeByIndex((*bead)->getResnr())->AddBead((*bead)->getId(), (*bead)->getName());

        MoleculeByIndex((*bead)->getResnr())->AddBead((*bead), string("1:TRI:") + (*bead)->getName());
    }
    
    /// \todo sort beads in molecules that all beads are stored in the same order. This is needed for the mapping!
}

void Topology::CreateOneBigMolecule(string name)
{
    Molecule *mi = CreateMolecule(name);
    
    BeadContainer::iterator bead;    
    
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
        stringstream n("");
        n << (*bead)->getResnr() +1 << ":" <<  _residues[(*bead)->getResnr()]->getName() << ":" << (*bead)->getName();
        //cout << n.str() << endl;
        mi->AddBead((*bead), n.str());
    }    
}

void Topology::Add(Topology *top)
{
    BeadContainer::iterator bead;
    ResidueContainer::iterator res;
    MoleculeContainer::iterator mol;
    
    int res0=ResidueCount();
    
    for(bead=top->_beads.begin(); bead!=top->_beads.end(); ++bead) {
        Bead *bi = *bead;
        BeadType *type =  GetOrCreateBeadType(bi->getType()->getName());
        CreateBead(bi->getSymmetry(), bi->getName(), type, bi->getResnr()+res0, bi->getMass(), bi->getQ());
    }
    
    for(res=top->_residues.begin();res!=top->_residues.end(); ++res) {
        CreateResidue((*res)->getName());
    }
  
    // \todo beadnames in molecules!!
    for(mol=top->_molecules.begin();mol!=top->_molecules.end(); ++mol) {
        Molecule *mi = CreateMolecule((*mol)->getName());
        for(int i=0; i<mi->BeadCount(); i++) {
            mi->AddBead(mi->getBead(i), "invalid");
        }
    }
}

void Topology::CopyTopologyData(Topology *top)
{
    BeadContainer::iterator it_bead;
    ResidueContainer::iterator it_res;
    MoleculeContainer::iterator it_mol;

    _bc->setBox(top->getBox());
    _time = top->_time;
    _step = top->_step;

    // cleanup old data
    Cleanup();

    // copy all residues
    for(it_res=top->_residues.begin();it_res!=top->_residues.end(); ++it_res) {
        CreateResidue((*it_res)->getName());
    }

    // create all beads
    for(it_bead=top->_beads.begin(); it_bead!=top->_beads.end(); ++it_bead) {
        Bead *bi = *it_bead;
        BeadType *type =  GetOrCreateBeadType(bi->getType()->getName());
        Bead *bn = CreateBead(bi->getSymmetry(), bi->getName(), type, bi->getResnr(), bi->getMass(), bi->getQ());
        bn->setOptions(bi->Options());
    }

    // copy all molecules
    for(it_mol=top->_molecules.begin();it_mol!=top->_molecules.end(); ++it_mol) {
        Molecule *mi = CreateMolecule((*it_mol)->getName());
        for(int i=0; i<(*it_mol)->BeadCount(); i++) {
            int beadid = (*it_mol)->getBead(i)->getId();
            mi->AddBead(_beads[beadid], (*it_mol)->getBeadName(i));
        }
    }
    // TODO: copy interactions
    //InteractionContainer::iterator it_ia;
    //for(it_ia=top->_interaction.begin();it_ia=top->_interactions.end();++it_ia) {

    //}
}


void Topology::RenameMolecules(string range, string name)
{
    RangeParser rp;
    RangeParser::iterator i;
    
    rp.Parse(range);
    for(i=rp.begin();i!=rp.end();++i) {
        if((unsigned int)*i > _molecules.size())
            throw runtime_error(string("RenameMolecules: num molecules smaller than"));
        getMolecule(*i-1)->setName(name);
    }
}

void Topology::RenameBeadType(string name, string newname)
{
    BeadContainer::iterator bead;
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
      BeadType *type =  GetOrCreateBeadType((*bead)->getType()->getName());
      if (wildcmp(name.c_str(),(*bead)->getType()->getName().c_str())) {
	type->setName(newname);
      }
    }
}

void Topology::SetBeadTypeMass(string name, double value)
{
    BeadContainer::iterator bead;
    for(bead=_beads.begin(); bead!=_beads.end(); ++bead) {
      if (wildcmp(name.c_str(),(*bead)->getType()->getName().c_str())) {
	(*bead)->setMass(value);
      }
    }

}

void Topology::CheckMoleculeNaming(void)
{
    map<string,int> nbeads;

    for(MoleculeContainer::iterator iter = _molecules.begin(); iter!=_molecules.end(); ++iter) {
        map<string,int>::iterator entry = nbeads.find((*iter)->getName());
        if(entry != nbeads.end()) {
            if(entry->second != (*iter)->BeadCount())
                throw runtime_error("There are molecules which have the same name but different number of bead "
                        "please check the section manual topology handling in the votca manual");
            continue;
        }
        nbeads[(*iter)->getName()] = (*iter)->BeadCount();
    }
}


void Topology::AddBondedInteraction(Interaction *ic)
{
    map<string,int>::iterator iter;
    iter = _interaction_groups.find(ic->getGroup());
    if(iter!=_interaction_groups.end())
        ic->setGroupId((*iter).second);
    else {
        int i= _interaction_groups.size();
        _interaction_groups[ic->getGroup()] = i;
        ic->setGroupId(i);
    }
    _interactions.push_back(ic);
    _interactions_by_group[ic->getGroup()].push_back(ic);
}

std::list<Interaction *> Topology::InteractionsInGroup(const string &group)
{
    map<string, list<Interaction*> >::iterator iter;
    iter = _interactions_by_group.find(group);
    if(iter == _interactions_by_group.end())
        return list<Interaction *>();
    return iter->second;
}


BeadType *Topology::GetOrCreateBeadType(string name)
{
    map<string, int>::iterator iter;
    
    iter = _beadtype_map.find(name);
    if(iter == _beadtype_map.end()) {
        BeadType *bt = new BeadType(this, _beadtypes.size(), name);
        _beadtypes.push_back(bt);
        _beadtype_map[name] = bt->getId();
        return bt;
    }
    else {
        return _beadtypes[(*iter).second];
    }
    return NULL;
}

vec Topology::BCShortestConnection(const vec &r_i, const vec &r_j) const
{
    return _bc->BCShortestConnection(r_i, r_j);
}

vec Topology::getDist(int bead1, int bead2) const
{
    return BCShortestConnection(
            getBead(bead1)->getPos(),
            getBead(bead2)->getPos());
}

double Topology::BoxVolume()
{
    return _bc->BoxVolume();
}

void Topology::RebuildExclusions()
{
    _exclusions.CreateExclusions(this);
}

BoundaryCondition::eBoxtype Topology::autoDetectBoxType(const matrix &box) {
    // set the box type to OpenBox in case "box" is the zero matrix,
    // to OrthorhombicBox in case "box" is a diagonal matrix,
    // or to TriclinicBox otherwise
    if(box.get(0,0)==0 && box.get(0,1)==0 && box.get(0,2)==0 &&
       box.get(1,0)==0 && box.get(1,1)==0 && box.get(1,2)==0 &&
       box.get(2,0)==0 && box.get(2,1)==0 && box.get(2,2)==0) {
        //cout << "box open\n";
        return BoundaryCondition::typeOpen;
    }
    else
    if(box.get(0,1)==0 && box.get(0,2)==0 &&
       box.get(1,0)==0 && box.get(1,2)==0 &&
       box.get(2,0)==0 && box.get(2,1)==0) {
        //cout << "box orth\n";
        return BoundaryCondition::typeOrthorhombic;
    }
    else {
        //cout << "box tric\n";
        return BoundaryCondition::typeTriclinic;
    }
    return BoundaryCondition::typeOpen;
}

double Topology::ShortestBoxSize()
{
    vec _box_a = getBox().getCol(0);
    vec _box_b = getBox().getCol(1);
    vec _box_c = getBox().getCol(2);

    // create plane normals
    vec _norm_a = _box_b ^ _box_c;
    vec _norm_b = _box_c ^ _box_a;
    vec _norm_c = _box_a ^ _box_b;

    _norm_a.normalize();
    _norm_b.normalize();
    _norm_c.normalize();

    double la = _box_a * _norm_a;
    double lb = _box_b * _norm_b;
    double lc = _box_c * _norm_c;

    return min(la, min(lb, lc));
}

}}
