#include "jcalc.h"


void JCalc::Init(string filename){

    xmlDocPtr doc;
    xmlNodePtr node;
    xmlChar *key;

    // open the xml file
    //cout << filename << endl;
    doc = xmlParseFile(filename.c_str());
    if(doc == NULL)
        throw "Error on open crgunittype list: " + filename;

    node = xmlDocGetRootElement(doc);

    if(node == NULL) {
        xmlFreeDoc(doc);
        throw "Error, empty xml document: " + filename;
    }

    if(xmlStrcmp(node->name, (const xmlChar *)"crgunit_type")) {
        xmlFreeDoc(doc);
        throw "Error, xml file not labeled crgunit_type: " + filename;
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

