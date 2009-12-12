#include "jcalc.h"


void JCalc::Init(string filename){

    xmlDocPtr doc;
    xmlNodePtr node;
    xmlChar *key;

    // open the xml file
    //cout << filename << endl;
    doc = xmlParseFile(filename.c_str());
    if(doc == NULL)
        throw runtime_error("Error on open crgunittype list: " + filename);

    node = xmlDocGetRootElement(doc);

    if(node == NULL) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, empty xml document: " + filename);
    }

    if(xmlStrcmp(node->name, (const xmlChar *)"crgunit_type")) {
        xmlFreeDoc(doc);
        throw runtime_error("Error, xml file not labeled crgunit_type: " + filename);
    }
    // parse xml tree
    for(node = node->xmlChildrenNode; node != NULL; node = node->next) {
            key = xmlNodeListGetString(doc, node->xmlChildrenNode, 1);
            if (!xmlStrcmp( node->name, (const xmlChar *) "ChargeUnitType")){
                ParseCrgUnitType(doc, node->xmlChildrenNode);
            }

    }
    #ifdef DEBUG
    cout << "Finished reading the list of charges, we got  " << _listCrgUnitType.size()
         << "types of charges " <<endl;
    #endif
    // initialise the list of rate calculators
//    init_listfock();

}


void JCalc::ParseCrgUnitType(xmlDocPtr doc, xmlNodePtr cur ){
    xmlChar *key;

    char * namecoord;
    char * nameorb;
    char * namecrg;
    char * nameneutr;
    double reorg;
    double energy;
    // int transorb;
    vector <unsigned int> transorbs;
    string name;
    string beadconj;
    string molname;

    vector < vector <int> > list_atoms_monomer;
    vector < vector <double> > list_weights_monomer;

    for(; cur != NULL; cur = cur->next) {
            if (!xmlStrcmp(cur->name,(const xmlChar *) "posname")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                namecoord = reinterpret_cast<char *>(key);
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "orbname")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                nameorb = reinterpret_cast<char *>(key);
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "nameneutr")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                nameneutr = reinterpret_cast<char *>(key);
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "namecrg")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                namecrg = reinterpret_cast<char *>(key);
            }
            /*else if (!xmlStrcmp(cur->name, (const xmlChar *) "transorb")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                transorb = atoi ( reinterpret_cast<char *>(key));
            }*/
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "transorb")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string s( reinterpret_cast<char *>(key));
                Tokenizer tok(s, " \n\t");
                Tokenizer::iterator beg=tok.begin();
                unsigned int value;
                while(beg!=tok.end()){
                    value = atoi((*beg).c_str());
                    transorbs.push_back(value);
                    ++beg;
                }
                clog << "Orbitals relevant for charge transport: ";
                for(int i=0; i<transorbs.size(); i++){
                    clog << transorbs[i] << " ";
                }
                clog << endl;
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "reorg")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                reorg = atof ( reinterpret_cast<char *>(key) );
            }
            else if (!xmlStrcmp(cur->name,  (const xmlChar *) "beadconj")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                beadconj = reinterpret_cast<char *>(key);
            }
            else if (!xmlStrcmp(cur->name,  (const xmlChar *) "name")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                name = reinterpret_cast<char *>(key);
            }
            else if (!xmlStrcmp(cur->name,  (const xmlChar *) "molname")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                molname = reinterpret_cast<char *>(key);
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "energy")){
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                energy = atof ( reinterpret_cast<char *>(key) ) *27.21;
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "monomer_atom_map")){
                #ifdef DEBUG
                cout << "readidng the monomer map" <<endl;
                #endif
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string all_monomer (reinterpret_cast<char *> (key));

                for(string::iterator c = all_monomer.begin(); c!=all_monomer.end(); ++c){
                    if(*c == '\n') *c = ' ';
                    if(*c == '\t') *c = ' ';
                }

                Tokenizer tok1(all_monomer, ":");
                Tokenizer::iterator it_mon;
                for (it_mon= tok1.begin(); it_mon != tok1.end() ; ++it_mon){
                    Tokenizer tok2 (*it_mon, " ");
                    Tokenizer::iterator it_at;
                    vector <int> list_atoms;
                    for (it_at = tok2.begin() ; it_at != tok2.end(); ++it_at){
                        list_atoms.push_back (lexical_cast<int> (*it_at));
                    }
                    list_atoms_monomer.push_back(list_atoms);
                }
            }
            else if (!xmlStrcmp(cur->name, (const xmlChar *) "monomer_atom_weights")){
                #ifdef DEBUG
                cout << "readidng the monomer map" <<endl;
                #endif
                key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                string all_monomer (reinterpret_cast<char *> (key));
                for(string::iterator c = all_monomer.begin(); c!=all_monomer.end(); ++c){
                    if(*c == '\n') *c = ' ';
                    if(*c == '\t') *c = ' ';
                }
                Tokenizer tok1(all_monomer, ":");
                Tokenizer::iterator it_mon;
                for (it_mon= tok1.begin(); it_mon != tok1.end() ; ++it_mon){
                    Tokenizer tok2 (*it_mon, " ");
                    Tokenizer::iterator it_at;
                    vector <double> list_weights;
                    for (it_at = tok2.begin() ; it_at != tok2.end(); ++it_at){
                        list_weights.push_back (lexical_cast<double> (*it_at));
                    }
                    list_weights_monomer.push_back(list_weights);
                }
            }
    }
    CrgUnitType* crgunittype = new CrgUnitType(namecoord, nameorb, nameneutr, namecrg, reorg, energy, transorbs, _listCrgUnitType.size() ,
            molname, name, list_atoms_monomer, list_weights_monomer);
    _mapCrgUnitByName.insert( make_pair(name, crgunittype));

    clearListList(list_atoms_monomer);
    _listCrgUnitType.push_back(crgunittype);
}

JCalc::JCalcData * JCalc::InitJCalcData(CrgUnitType * mol1, CrgUnitType *mol2 ){
    JCalcData * data = new JCalcData;
    
    data->_type1 = mol1;
    data->_type2 = mol2;

    data->_orblabels.first = data->_type1 -> GetTransOrbs();
    data->_orblabels.second = data->_type2 -> GetTransOrbs();

    int nrorbs = data->_orblabels.first.size();
    for(int i=0; i<nrorbs; i++){
        if (data->_orblabels.first[i] != i && data->_orblabels.second[i] != i){
            cout << "orblabels: " << data->_orblabels.first[i] << " " << data->_orblabels.second[i] << endl;
             throw "Error in RateCalculator, the charge unit types do not have stripped orbitals";
        }
    }
    //initialise the first copy of the molecules + orbitals
    data->_mol1.define_bs(data->_indo);
    data->_mol1.cp_atompos(mol1->GetCrgUnit() );
    data->_mol1.cp_atoms  (mol1->GetCrgUnit() );
    data->_orb1.init_orbitals_stripped(mol1->GetOrb(), nrorbs);
    data->_mol1.assign_orb(&data->_orb1);
    data->_mol1.cp_crg(mol1->GetCrgUnit());

 /***   cout << "crgunit type: " <<endl;
    cout << mol1->GetCrgUnit().getN() <<endl;
    cout << "data type: " <<endl;
    cout << data->_mol1.getN() <<endl;**/
    //inititalise the second copy of molecules + orbitals
    data->_mol2.define_bs(data->_indo);
    data->_mol2.cp_atompos(mol2->GetCrgUnit() );
    data->_mol2.cp_atoms  (mol2->GetCrgUnit() );
    data->_orb2.init_orbitals_stripped(mol2->GetOrb(), nrorbs);
    data->_mol2.assign_orb(&data->_orb2);
    data->_mol2.cp_crg(mol2->GetCrgUnit());

    /***cout << "crgunit type: " <<endl;
    cout << mol2->GetCrgUnit().getN() <<endl;
    cout << "data type: " <<endl;
    cout << data->_mol2.getN() <<endl;
    **///
    // we have stripped the orbs to the bone
    for(int i=0; i < data->_orblabels.first.size(); i++){
        data->_orblabels.first[i]  = i;
        data->_orblabels.second[i] = i;
    }

    //initialise the fock matrix
    data->_fock.init(data->_mol1, data->_mol2);

    _maplistfock.insert(make_pair(make_pair(mol1, mol2) , data ));
    return data;
}

vector <double> JCalc::GetJ (CrgUnit & one, CrgUnit & two) {
    if ( one.GetTypeID() > two.GetTypeID()){
        return GetJ(two, one);
    }

    JCalcData * jdata;
    map <pair<CrgUnitType *, CrgUnitType *> , JCalcData *>::iterator itm=
             _maplistfock.find(make_pair(one.GetType(), two.GetType()));
     if (itm == _maplistfock.end() ){
         jdata = InitJCalcData(one.GetType(), two.GetType());
///         cout << "after creating" <<endl;
///         cout<< jdata->_mol1.getN()<<endl;
     }
     else{
         jdata = itm->second;
     }

    //calculate J
    #ifdef DEBUG
    cout << one.GetId() << "is the first molecule I am rotatin and the second is: " <<  two.GetId() <<endl;
    #endif

    one.rot_two_mol(two, jdata->_mol1, jdata->_mol2);

    // calculating orbitals for all pairs

    vector < double > Js; // vector with all Js
    vector < pair <int,int> > _pois; // vector with all pairs of interest
    pair <int, int> poi; // pair of interest
    for(int i=0; i<jdata->_orblabels.first.size(); i++){
        for(int j=0; j<jdata->_orblabels.first.size(); j++){
            poi.first = jdata->_orblabels.first[i];
            poi.second = jdata->_orblabels.second[j];
            Js.push_back(jdata->_fock.calcJ(poi));
            _pois.push_back(poi);
        }
    }

    //copy molecules back
    jdata->_mol1.cp_atompos(jdata->_type1->GetCrgUnit());
    jdata->_mol2.cp_atompos(jdata->_type2->GetCrgUnit());

    jdata->_mol1.cp_orb(jdata->_type1->GetCrgUnit(), jdata->_orblabels.first);
    jdata->_mol2.cp_orb(jdata->_type2->GetCrgUnit(), jdata->_orblabels.second);

    return Js;
}

CrgUnitType * JCalc::GetCrgUnitTypeByName(string name){
    map <string, CrgUnitType *>::iterator ittype = this->_mapCrgUnitByName.find(name);
    if (ittype == _mapCrgUnitByName.end()){
        throw runtime_error("Cannot find Crg unit type with name" + name);
    }
    CrgUnitType *type = ittype->second;
    return type;
}

CrgUnit JCalc::DefineCrgUnit(vec pos, matrix orient, string name){
    CrgUnitType *type = GetCrgUnitTypeByName(name);

    vec plane1 = orient.getCol(1);
    vec norm1 = orient.getCol(2);
    return CrgUnit::CrgUnit(pos, plane1, norm1, type);
}