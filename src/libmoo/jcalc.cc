#include "jcalc.h"

JCalc::~JCalc()
{
    safe_delete(_listCrgUnitType);
    _mapCrgUnitByName.clear();

    map <pair<CrgUnitType *, CrgUnitType *>, JCalcData *>::iterator itm = _maplistfock.begin();
    for (; itm != _maplistfock.end(); itm++)
        delete itm->second;
    _maplistfock.clear();
}

/// TODO: rewrite this using Property class
/*void JCalc::Initialize(const string &filename)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    xmlChar *key;

    // open the xml file
    //cout << filename << endl;
    doc = xmlParseFile(filename.c_str());
    if (doc == NULL)
        throw runtime_error("Error on open crgunittype list: " + filename);

    node = xmlDocGetRootElement(doc);

    if (node == NULL) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, empty xml document: " + filename);
    }

    if (xmlStrcmp(node->name, (const xmlChar *) "crgunit_type")) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, xml file not labeled crgunit_type: " + filename);
    }
    // parse xml tree
    for (node = node->xmlChildrenNode; node != NULL; node = node->next) {
        key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
        if (!xmlStrcmp(node->name, (const xmlChar *) "ChargeUnitType")) {
            ParseCrgUnitType(doc, node->xmlChildrenNode);
        }
    }
}*/

void JCalc::Initialize(const string& filename){
    load_property_from_xml(_options, filename.c_str());

    ParseCrgUnitTypes(_options.get("crgunit_type"));
}

void JCalc::ParseCrgUnitTypes(Property &options){
    string namecoord;
    string nameorb;
    string namecrg;
    string nameneutr;
    double reorg;
    double energy;
    vector <int> transorbs;
    string beadconj;
    string name;
    string molname;
    string namebasis;
    vector < vector <int> > list_atoms_monomer;
    vector < vector <double> > list_weights_monomer;

    list<Property *> types = options.Select("ChargeUnitType");
    list<Property *>::iterator iter;

    for (iter = types.begin();iter != types.end(); ++iter){
        namecoord = (*iter)->get("posname").as<string>();
        nameorb = (*iter)->get("orbname").as<string>();
        bool estatics = (*iter)->exists("nameneutr");
        if(estatics==true){
            nameneutr = (*iter)->get("nameneutr").as<string>();
            namecrg = (*iter)->get("namecrg").as<string>();
        }
        reorg = (*iter)->get("reorg").as<double>();
        energy = (*iter)->get("energy").as<double>();
        transorbs = (*iter)->get("transorb").as<vector <int> >();
        beadconj = (*iter)->get("beadconj").as<string>();
        name = (*iter)->get("name").as<string>();
        molname = (*iter)->get("molname").as<string>();
        namebasis = (*iter)->get("basisset").as<string>();

        string all_monomer = (*iter)->get("monomer_atom_map").as<string>();
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
            list_atoms_monomer.push_back(list_atoms);
        }

        string all_weights = (*iter)->get("monomer_atom_weights").as<string>();
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
    }
    CrgUnitType* crgunittype = new CrgUnitType(namecoord.c_str(), nameorb.c_str(),
            nameneutr.c_str(), namecrg.c_str(), namebasis,
            reorg, energy, transorbs, _listCrgUnitType.size(),
            molname, name, list_atoms_monomer, list_weights_monomer);
    _mapCrgUnitByName.insert(make_pair(name, crgunittype));

    clearListList(list_atoms_monomer);
    clearListList(list_weights_monomer);
    _listCrgUnitType.push_back(crgunittype);
}

/// TODO: rewrite this using Property class

/*void JCalc::ParseCrgUnitType(xmlDocPtr doc, xmlNodePtr cur)
{
    xmlChar *key;

    string namecoord;
    string nameorb;
    string namecrg;
    string nameneutr;
    double reorg;
    double energy;
    // int transorb;
    vector < int> transorbs;
    string name;
    string beadconj;
    string molname;
    string basisset= string("INDO");

    vector < vector <int> > list_atoms_monomer;
    vector < vector <double> > list_weights_monomer;

    for (; cur != NULL; cur = cur->next) {
        if (!xmlStrcmp(cur->name, (const xmlChar *) "posname")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            namecoord = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "orbname")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            nameorb = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "nameneutr")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            nameneutr = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "namecrg")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            namecrg = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "transorb")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string s(reinterpret_cast<char *> (key));
            Tokenizer tok(s, " \n\t");
            Tokenizer::iterator beg = tok.begin();
            int value;
            while (beg != tok.end()) {
                value = atoi((*beg).c_str());
                transorbs.push_back(value);
                ++beg;
            }
            clog << "Orbitals relevant for charge transport: ";
            for (int i = 0; i < transorbs.size(); i++) {
                clog << transorbs[i] << " ";
            }
            clog << endl;
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "reorg")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            reorg = atof(reinterpret_cast<char *> (key));
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "beadconj")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            beadconj = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "name")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            name = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "molname")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            molname = reinterpret_cast<char *> (key);
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "energy")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            energy = atof(reinterpret_cast<char *> (key)) *27.21;
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "basisset")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            basisset =string (reinterpret_cast<char *> (key)) ;
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "monomer_atom_map")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string all_monomer(reinterpret_cast<char *> (key));

            for (string::iterator c = all_monomer.begin(); c != all_monomer.end(); ++c) {
                if (*c == '\n') *c = ' ';
                if (*c == '\t') *c = ' ';
            }

            Tokenizer tok1(all_monomer, ":");
            Tokenizer::iterator it_mon;
            for (it_mon = tok1.begin(); it_mon != tok1.end(); ++it_mon) {
                Tokenizer tok2(*it_mon, " ");
                Tokenizer::iterator it_at;
                vector <int> list_atoms;
                for (it_at = tok2.begin(); it_at != tok2.end(); ++it_at) {
                    list_atoms.push_back(lexical_cast<int> (*it_at));
                }
                list_atoms_monomer.push_back(list_atoms);
            }
        }
        else if (!xmlStrcmp(cur->name, (const xmlChar *) "monomer_atom_weights")) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string all_monomer(reinterpret_cast<char *> (key));
            for (string::iterator c = all_monomer.begin(); c != all_monomer.end(); ++c) {
                if (*c == '\n') *c = ' ';
                if (*c == '\t') *c = ' ';
            }
            Tokenizer tok1(all_monomer, ":");
            Tokenizer::iterator it_mon;
            for (it_mon = tok1.begin(); it_mon != tok1.end(); ++it_mon) {
                Tokenizer tok2(*it_mon, " ");
                Tokenizer::iterator it_at;
                vector <double> list_weights;
                for (it_at = tok2.begin(); it_at != tok2.end(); ++it_at) {
                    list_weights.push_back(lexical_cast<double> (*it_at));
                }
                list_weights_monomer.push_back(list_weights);
            }
        }
    }
    CrgUnitType* crgunittype = new CrgUnitType(namecoord.c_str(), nameorb.c_str(), 
            nameneutr.c_str(), namecrg.c_str(), basisset,
            reorg, energy, transorbs, _listCrgUnitType.size(),
            molname, name, list_atoms_monomer, list_weights_monomer);
    _mapCrgUnitByName.insert(make_pair(name, crgunittype));

    clearListList(list_atoms_monomer);
    _listCrgUnitType.push_back(crgunittype);
}*/

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
    for (int i = 0; i < data->_orblabels.first.size(); i++) {
        data->_orblabels.first[i] = i;
    }
    for (int i =0;i < data->_orblabels.second.size(); i++) {
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
    if (one.getType()->getId() > two.getType()->getId())
        return CalcJ(two, one);
    
    JCalcData * jdata = getJCalcData(one, two);

    //calculate J
    one.rot_two_mol(two, jdata->_mol1, jdata->_mol2);

    // calculating orbitals for all pairs

    vector < double > Js; // vector with all Js
    vector < pair <int, int> > _pois; // vector with all pairs of interest
    pair <int, int> poi; // pair of interest
    for (int i = 0; i < jdata->_orblabels.first.size(); i++) {
        for (int j = 0; j < jdata->_orblabels.first.size(); j++) {
            poi.first = jdata->_orblabels.first[i];
            poi.second = jdata->_orblabels.second[j];
            Js.push_back(jdata->_fock.calcJ(poi));
            _pois.push_back(poi);
        }
    }
    return Js;
}

int JCalc::WriteProJ(CrgUnit & one, CrgUnit & two)
{
    // only write them when one.id < two.id
    if (one.getType()->getId() < two.getType()->getId())
        return 1;

    //rotate the molecule and orbitals
    JCalcData * jdata = getJCalcData(one, two);
    one.rot_two_mol(two, jdata->_mol1, jdata->_mol2);


    ofstream out;
    string name1, nameorb1, name2, nameorb2, namedim;
    name1 = lexical_cast<string>(one.getId()) + ".xyz";
    name2 = lexical_cast<string>(two.getId()) + ".xyz";
    nameorb1 = lexical_cast<string>(one.getId()) + ".fort.7";
    nameorb2 = lexical_cast<string>(two.getId()) + ".fort.7";
    namedim = lexical_cast<string>(one.getId())+"and"+
            lexical_cast<string>(two.getId()) +".com";

    
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
           "#p b3lyp/TZVP-6D guess(cards) scf(MaxCycle=1, Conver=1) "<<
            "Nosymm IOp(3/33=1) IOp(5/32=1) punch=mo" << '\n' <<
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

double JCalc::EstaticDifference(CrgUnit & crged, CrgUnit & neutr)
{

    // we need to calculate E:CRG-NEUTR - E:NEUTR-NEUTR
    // the difficulty is that ratecalculator
    // exists only for type1 < type2
    // this requires to check the typeIDs
    JCalcData * jdata = getJCalcData(crged, neutr);

    if (crged.getType()->getId() > neutr.getType()->getId()) {

        crged.rot_two_mol(neutr, jdata->_mol2, jdata->_mol1);

        double nrg = jdata->_mol2.V_nrg_crg_neutr(jdata->_mol1) -
                jdata->_mol2.V_nrg_neutr_neutr(jdata->_mol1);
        //copy molecules back
        jdata->_mol1.cp_atompos(jdata->_type1->GetCrgUnit());
        jdata->_mol2.cp_atompos(jdata->_type2->GetCrgUnit());
        jdata->_mol1.cp_orb(jdata->_type1->GetCrgUnit(), jdata->_orblabels.first);
        jdata->_mol2.cp_orb(jdata->_type2->GetCrgUnit(), jdata->_orblabels.second);
        return unit<hartree,eV>::to(nrg);
    }
    else {
        crged.rot_two_mol(neutr, jdata->_mol1, jdata->_mol2);

        double nrg = jdata->_mol1.V_nrg_crg_neutr(jdata->_mol2) - jdata->_mol1.V_nrg_neutr_neutr(jdata->_mol2);
        //copy molecules back
        jdata->_mol1.cp_atompos(jdata->_type1->GetCrgUnit());
        jdata->_mol2.cp_atompos(jdata->_type2->GetCrgUnit());
        jdata->_mol1.cp_orb(jdata->_type1->GetCrgUnit(), jdata->_orblabels.first);
        jdata->_mol2.cp_orb(jdata->_type2->GetCrgUnit(), jdata->_orblabels.second);
        return unit<hartree,eV>::to(nrg);
    }
}

CrgUnitType * JCalc::GetCrgUnitTypeByName(string name)
{
    map <string, CrgUnitType *>::iterator ittype = this->_mapCrgUnitByName.find(name);
    if (ittype == _mapCrgUnitByName.end()) {
        throw runtime_error("Cannot find Crg unit type with name" + name);
    }
    CrgUnitType *type = ittype->second;
    return type;
}
