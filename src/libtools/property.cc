/* 
 * File:   property.cc
 * Author: ruehle
 *
 * Created on November 13, 2008, 5:54 PM
 */

#include <libxml/parser.h>
#include <iostream>
#include "property.h"


Property &Property::get(const string &key)
{
    Tokenizer tok(key, ".");
    Tokenizer::iterator n;
    
    n = tok.begin();
    
    map<string, Property>::iterator iter;
    iter = _map.find(*n);
    if(iter == _map.end())
        throw string("property not found: ") + key;
    
    Property *p(&((*iter).second));
    ++n;
    try {
        for(; n!=tok.end(); ++n)
        p = &p->get(*n);
    }
    catch(string err) {
        throw string("property not found: ") + key;
    }
    
    return *p;
}

void Property::PrintNode(std::ostream &out, const string &prefix, Property &p)
{
    
    map<string, Property>::iterator iter;
    if((p._value != "") || p.HasChilds())
        out << prefix << " = " << p._value << endl;
    for(iter = p._map.begin(); iter!=p._map.end(); ++iter) {
        if(prefix=="") 
            PrintNode(out, prefix + (*iter).first, (*iter).second);
        else
            PrintNode(out, prefix + "." + (*iter).first, (*iter).second);
    }
}

static void parse_node(Property &p, xmlDocPtr doc, xmlNodePtr cur)
{
    xmlChar *key;
    
    key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
    
    string value("");
    if(key) value = (const char *)key;
    
    Property &np = p.set((const char *)cur->name, value);
    xmlFree(key);    

    for(xmlNodePtr node = cur->xmlChildrenNode; node != NULL; node = node->next) {
        if (node->type == XML_ELEMENT_NODE)
            parse_node(np, doc, node);                    
    }
}

bool load_property_from_xml(Property &p, string filename)
{
    xmlDocPtr doc;
    xmlNodePtr node;
    xmlChar *key;
    
    doc = xmlParseFile(filename.c_str());
    if(doc == NULL) 
        throw std::ios_base::failure("Error on open xml bead map: " + filename);
    
    node = xmlDocGetRootElement(doc);
    
    if(node == NULL) {
        xmlFreeDoc(doc);
        throw std::invalid_argument("Error, empty xml document: " + filename);
    }
    
    //if(xmlStrcmp(node->name, (const xmlChar *)"cg_molecule")) {
    //    xmlFreeDoc(doc);
    //    throw std::invalid_argument("Error, in xml file: " + filename);
    //}
    
    // parse xml tree
    parse_node(p, doc, node);
    
    // free the document
    xmlFreeDoc(doc);
    xmlCleanupParser();
}

std::ostream &operator<<(std::ostream &out, Property& p)
{
      Property::PrintNode(out, "", p);
      return out;
}
