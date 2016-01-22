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

#include <iostream>
#include <votca/xtp/jcalc.h>
#include <votca/tools/globals.h>

namespace votca { namespace xtp {

JCalc::~JCalc()
{
    safe_delete(_listCrgUnitType);
    _mapCrgUnitByName.clear();

    map <pair<CrgUnitType *, CrgUnitType *>, JCalcData *>::iterator itm = _maplistfock.begin();
    for (; itm != _maplistfock.end(); itm++)
        delete itm->second;
    _maplistfock.clear();
}

void JCalc::Initialize(const string& filename){
    load_property_from_xml(_options, filename.c_str());

    ParseCrgUnitTypes(_options.get("segments"));
}

void JCalc::ParseCrgUnitTypes(Property &options){
    string namecoord;
    string nameorb;
    string namecrg;
    string nameneutr;
    vector <int> transorbs;
   
    string name;
    string molname;
    string namebasis;
    vector < vector <int> > list_atoms_monomer;
    vector < vector <double> > list_weights_monomer;

    list<Property *> types = options.Select("segment");
    list<Property *>::iterator iter;

    for (iter = types.begin();iter != types.end(); ++iter){
        namecoord = (*iter)->get("coordinates").as<string>();

	bool estatics = (*iter)->exists("qneutral");
        if(estatics==true){
            nameneutr = (*iter)->get("qneutral").as<string>();
            namecrg = (*iter)->get("qcharged").as<string>();
        }
        
        if((*iter)->exists("torbital")) transorbs = (*iter)->get("torbital").as<vector <int> >();
        for(unsigned int i=0; i<transorbs.size(); ++i)
            transorbs[i]--;
        
        if(transorbs.size() > 0) nameorb = (*iter)->get("orbitals").as<string>();

        name = (*iter)->get("name").as<string>();
        
        namebasis = "INDO";
        if((*iter)->exists("basisset") > 0) namebasis = (*iter)->get("basisset").as<string>();

        string all_monomer = (*iter)->get("map").as<string>();
        for (string::iterator c = all_monomer.begin(); c != all_monomer.end(); ++c) {
            if (*c == '\n') *c = ' ';
            if (*c == '\t') *c = ' ';
        }
        Tokenizer tok1(all_monomer, ":");
        Tokenizer::iterator it_mon;
        for (it_mon = tok1.begin(); it_mon != tok1.end(); ++it_mon) {
            vector <int> list_atoms;
            Tokenizer tok2(*it_mon, " ");
            tok2.ConvertToVector<int>(list_atoms);
            for(unsigned int i=0; i<list_atoms.size(); ++i) {
                if(list_atoms[i] == 0)
                    throw std::runtime_error("An atom with index=0 was specified in listcharges. Counting starts at 1!");

                list_atoms[i]--;
            }
            list_atoms_monomer.push_back(list_atoms);
        }

        string all_weights = (*iter)->get("weights").as<string>();
        for (string::iterator c = all_weights.begin(); c != all_weights.end(); ++c) {
            if (*c == '\n') *c = ' ';
            if (*c == '\t') *c = ' ';
        }
        Tokenizer tok3(all_weights, ":");
        for (it_mon = tok3.begin(); it_mon != tok3.end(); ++it_mon) {
            vector <double> list_weights;
            Tokenizer tok2(*it_mon, " ");
            tok2.ConvertToVector<double>(list_weights);
            list_weights_monomer.push_back(list_weights);
        }

        CrgUnitType* crgunittype = new CrgUnitType(namecoord.c_str(), nameorb.c_str(),
                nameneutr.c_str(), namecrg.c_str(), namebasis,
                transorbs, _listCrgUnitType.size(),
                name, list_atoms_monomer, list_weights_monomer);
        crgunittype->setOptions(*iter);
        _mapCrgUnitByName.insert(make_pair(name, crgunittype));
        _listCrgUnitType.push_back(crgunittype);

        clearListList(list_atoms_monomer);
        clearListList(list_weights_monomer);
    }
}


JCalc::JCalcData * JCalc::InitJCalcData(CrgUnitType * mol1, CrgUnitType *mol2)
{
    JCalcData * data = new JCalcData;

    data->_type1 = mol1;
    data->_type2 = mol2;

    data->_orblabels.first = data->_type1 -> GetTransOrbs();
    data->_orblabels.second = data->_type2 -> GetTransOrbs();

    int nrorbs1 = data->_orblabels.first.size();
    int nrorbs2 = data->_orblabels.second.size();
    /*for (int i = 0; i < nrorbs; i++) {
        if (data->_orblabels.first[i] != i && data->_orblabels.second[i] != i) {
            cout << "orblabels: " << data->_orblabels.first[i] << " " << data->_orblabels.second[i] << endl;
            throw "Error in RateCalculator, the charge unit types do not have stripped orbitals";
        }
    }*/
    //initialise the first copy of the molecules + orbitals
    data->_mol1.define_bs(mol1->GetBS());
    data->_mol1.cp_atompos(mol1->GetCrgUnit());
    data->_mol1.cp_atoms(mol1->GetCrgUnit());
    data->_orb1.init_orbitals_stripped(mol1->GetOrb(), nrorbs1);
    data->_mol1.assign_orb(&data->_orb1);
    data->_mol1.cp_crg(mol1->GetCrgUnit());

    //inititalise the second copy of molecules + orbitals
    data->_mol2.define_bs(mol2->GetBS());
    data->_mol2.cp_atompos(mol2->GetCrgUnit());
    data->_mol2.cp_atoms(mol2->GetCrgUnit());
    data->_orb2.init_orbitals_stripped(mol2->GetOrb(), nrorbs2);
    data->_mol2.assign_orb(&data->_orb2);
    data->_mol2.cp_crg(mol2->GetCrgUnit());

    // we have stripped the orbs to the bone
    for (unsigned int i = 0; i < data->_orblabels.first.size(); i++) {
        data->_orblabels.first[i] = i;
    }
    for (unsigned int i =0;i < data->_orblabels.second.size(); i++) {
        data->_orblabels.second[i] = i;
    }

    //initialise the fock matrix
    data->_fock.init(data->_mol1, data->_mol2);

    // add a quick reference to find calculator for a pair
    _maplistfock.insert(make_pair(make_pair(mol1, mol2), data));

    return data;
}

JCalc::JCalcData * JCalc::getJCalcData(CrgUnit & one, CrgUnit & two)
{
    JCalcData * jdata;
    map <pair<CrgUnitType *, CrgUnitType *>, JCalcData *>::iterator itm =
            _maplistfock.find(make_pair(one.getType(), two.getType()));
    if (itm == _maplistfock.end()) {
        jdata = InitJCalcData(one.getType(), two.getType());
    }
    else {
        jdata = itm->second;
    }
    return jdata;
}

vector <double> JCalc::CalcJ(CrgUnit & one, CrgUnit & two)
{
    if(one.getType()->_transorbs.size() == 0)
        throw std::runtime_error("no orbital information found for conjugated segment " + one.getName());
    if(two.getType()->_transorbs.size() == 0)
        throw std::runtime_error("no orbital information found for conjugated segment " + two.getName());

    if (one.getType()->getId() > two.getType()->getId())
        return CalcJ(two, one);
    
    JCalcData * jdata = getJCalcData(one, two);

    //calculate J
    one.rot_two_mol(two, jdata->_mol1, jdata->_mol2);

    // calculating orbitals for all pairs

      if (tools::globals::verbose) {
           cout << "  id1:id2 [J] " << endl;
           cout << one.getType()->getId() << ":" 
                << two.getType()->getId() << " [";
       }    
    
    vector < double > Js; // vector with all Js
    vector < pair <int, int> > _pois; // vector with all pairs of interest
    pair <int, int> poi; // pair of interest
    for (unsigned int i = 0; i < jdata->_orblabels.first.size(); i++) {
        for (unsigned int j = 0; j < jdata->_orblabels.first.size(); j++) {
            poi.first = jdata->_orblabels.first[i];
            poi.second = jdata->_orblabels.second[j];
            Js.push_back(jdata->_fock.calcJ(poi));
            _pois.push_back(poi);
            if (tools::globals::verbose) { cout << jdata->_fock.calcJ(poi) << " " ; }
        }
        if (tools::globals::verbose) { cout << "]" << endl; }
    }
    
    
     
    
    int my = one.getType()->getId();
    cout << my << endl;
    
    return Js;
}

int JCalc::WriteProJ(CrgUnit & one, CrgUnit & two,string namedir)
{
    // only write them when one.id < two.id
    if (one.getType()->getId() < two.getType()->getId())
        return 1;

    //rotate the molecule and orbitals
    JCalcData * jdata = getJCalcData(one, two);
    one.rot_two_mol(two, jdata->_mol1, jdata->_mol2);


    ofstream out;
    string name1, nameorb1, name2, nameorb2, namedim;
    name1 =namedir + boost::lexical_cast<string>(one.getId()) + ".xyz";
    name2 = namedir +boost::lexical_cast<string>(two.getId()) + ".xyz";
    nameorb1 =namedir + boost::lexical_cast<string>(one.getId()) + ".fort.7";
    nameorb2 = namedir +boost::lexical_cast<string>(two.getId()) + ".fort.7";
    namedim = namedir +boost::lexical_cast<string>(one.getId())+"and"+
            boost::lexical_cast<string>(two.getId()) +".com";

    
    //write the info for molecule1
    out.open(name1.c_str());
    (jdata->_mol1).print(out);
    out.close();

    (jdata->_orb1).print_g03(nameorb1);
    //write the info for molecule2
    out.open(name2.c_str());
    (jdata->_mol2).print(out);
    out.close();

    (jdata->_orb2).print_g03(nameorb2);
    
    //write the input file for both
    orb dimerorb;
    dimerorb.dimerise_orbs((jdata->_orb1), (jdata->_orb2), (jdata->_mol1).n_el,(jdata->_mol2).n_el);

    out.open(namedim.c_str());
    // write header
    out << "%nproc=1" <<'\n' <<
           "%mem=1Gb" << '\n'<<
           "#p b3lyp/TZVP-6D guess(cards) scf(MaxCycle=1) "<<
            "Nosymm IOp(3/33=1) IOp(5/13=1) punch=mo" << '\n' <<
            '\n' <<
            "autogen" << '\n' << '\n'<<
            "0 1" <<'\n';
    (jdata->_mol1).print(out);
    (jdata->_mol2).print(out);
    //write blank line
    out << '\n';
    out.close();
    dimerorb.print_g03(namedim, string ("a"));
    return 0;
}

CrgUnitType * JCalc::GetCrgUnitTypeByName(string name)
{
    map <string, CrgUnitType *>::iterator ittype = this->_mapCrgUnitByName.find(name);
    if (ittype == _mapCrgUnitByName.end()) {
        throw runtime_error("Cannot find segment type with name: " + name);
    }
    CrgUnitType *type = ittype->second;
    return type;
}

}}
