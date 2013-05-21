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

#include <iostream>
#include <votca/tools/property.h>
#include <stdexcept>
#include <votca/tools/tokenizer.h>
#include <boost/algorithm/string.hpp>

#include <stdio.h>
#include <expat.h>
#include <string.h>
#include <fstream>
#include <string>
#include <stack>

namespace votca { namespace tools {

PropertyFormat XML(formatXML);
PropertyFormat TXT(formatTXT);
PropertyFormat T2T(formatT2T);
PropertyFormat LOG(formatLOG);
PropertyFormat TEX(formatTEX);

// ostream modifier defines the output format
const int Property::_format = std::ios_base::xalloc();  
    
Property &Property::get(const string &key)
{
    Tokenizer tok(key, ".");
    Tokenizer::iterator n;
    
    n = tok.begin();
    if(n==tok.end()) return *this;
    
    Property *p;
    map<string, Property*>::iterator iter;
    if(*n=="") {
        p = this;
    }
    else {
        iter = _map.find(*n);
    
        if(iter == _map.end())
            throw runtime_error("property not found: " + key);
    
        p = (((*iter).second));
    }
    ++n;
    try {
        for(; n!=tok.end(); ++n) {
            p = &p->get(*n);
        }
    }
    catch(string err) { // catch here to get full key in exception
        throw runtime_error("property not found: " + key);
    }
    
    return *p;
}

std::list<Property *> Property::Select(const string &filter)
{
    Tokenizer tok(filter, ".");
    
    std::list<Property *> selection;

    if(tok.begin()==tok.end()) return selection;
    
    selection.push_back(this);
        
    for (Tokenizer::iterator n = tok.begin();
            n != tok.end(); ++n) {
        std::list<Property *> childs;
        for (std::list<Property *>::iterator p = selection.begin();
                p != selection.end(); ++p) {
                for (list<Property>::iterator iter = (*p)->_properties.begin();
                    iter != (*p)->_properties.end(); ++iter) {
                    if (wildcmp((*n).c_str(), (*iter).name().c_str())) {
                        childs.push_back(&(*iter));
                    }
                }
        }
        selection = childs;        
    }

    return selection;
}

void Property::PrintNodeTXT(std::ostream &out, const string &prefix, Property &p)
{
    
    map<string, Property*>::iterator iter;
    if((p._value != "") || p.HasChilds()) {
        if ( prefix.find_first_not_of(' ') != std::string::npos )
            
        if((p._value).find_first_not_of("\t\n ") != std::string::npos)     
        out << prefix << " = ";
        
        if((p._value).find_first_not_of("\t\n ") != std::string::npos) 
            out << p._value << endl;
        
    }
        
    for(iter = p._map.begin(); iter!=p._map.end(); ++iter) {
        if(prefix=="") 
            PrintNodeTXT(out, prefix + (*iter).first, *(*iter).second);
        else
            PrintNodeTXT(out, prefix + "." + (*iter).first, *(*iter).second);
    }
}

void Property::PrintNodeXML(std::ostream &out, const string &prefix, Property &p, string offset)
{
    
    list<Property>::iterator iter;       
    map<string,string>::iterator ia;
    bool _endl = true;
    
    //cout << p._name << " Attributes: " << p._attributes.size() << endl;
    
    if( p.HasChilds() || p._value != "" ) {
        
        // the head node is always empty, do not print it
        if ( prefix.find_first_not_of(' ') != std::string::npos )  {
                out << offset << "<" << prefix ;
                // the attributes 
                for(ia = p._attributes.begin(); ia!=p._attributes.end(); ++ia) {
                        cout << " " << ia->first << "=\"" << ia->second << "\"" ;
                }
                out << ">";
        }
 
        if((p._value).find_first_not_of("\t\n ") != std::string::npos) {
            out << p._value;
            _endl = false;
        } else {
            out << endl;
            _endl = true;
        }
        
        for(iter = p._properties.begin(); iter!=p._properties.end(); ++iter) {
            offset += "\t";
            PrintNodeXML(out, (*iter)._name , (*iter), offset);
            offset.resize(offset.size()-1);           
        }
        
        if ( prefix.find_first_not_of(' ') != std::string::npos ) {
            if ( _endl ) {
                out << offset << "</" << prefix << ">"  << endl;
            } else {
                out << "</" << prefix << ">"  << endl;
            }
        }
    }
}

void Property::PrintNodeT2T(std::ostream &out, const string &prefix, Property &p) {
    out << "T2T format is not implemented";
}

void Property::PrintNodeLOG(std::ostream &out, const string &prefix, Property &p) {
    out << "LOG format is not implemented";
}

void Property::PrintNodeTEX(std::ostream &out, const string &prefix, Property &p, int offset) {

    list<Property>::iterator iter;       
     
    if((p._value != "") || p.HasChilds()) {
        
         if((p._value).find_first_not_of("\t\n ") != std::string::npos ) {
            string _tex_name = boost::replace_all_copy( p._name, "_", "\\_" );
            out << " \\hspace{" << offset << "pt} "
                << "\\hypertarget{" << prefix << "}"
                <<  "{" << _tex_name << "}" 
                << " & " <<  p._attributes["help"] << "\\\\" << endl;
         }
    }

        
    for(iter = p._properties.begin(); iter!=p._properties.end(); ++iter) {
        if(prefix=="") {
            offset += 10;
            PrintNodeTEX(out, prefix + (*iter)._name, (*iter), offset);
            offset -= 10;
        } else
            PrintNodeTEX(out, prefix + "." + (*iter)._name, (*iter), offset);
    }        
}

std::ostream &operator<<(std::ostream &out, Property& p)
{
    if (!out.good())
        return out;

    std::ostream::sentry sentry(out);

    if(sentry)
    {
        switch(out.iword(p.GetFormat()))
        {
        default:
            p.PrintNodeTXT(out, "", p);
        case formatXML:
            //cout << p._name << " " << p._value << " " << p._path << endl;
            p.PrintNodeXML(out, "", p, "");
            break;
        case formatTXT:
            p.PrintNodeTXT(out, "", p);
            break;
        case formatT2T:
            p.PrintNodeT2T(out, "", p);
            break;
        case formatLOG:
            p.PrintNodeLOG(out, "", p);
            break;
        case formatTEX:
            string name = p.begin()->_name;
            string label = (p.get(name))._attributes["label"];
            string help = (p.get(name))._attributes["help"];
            
            out << "\\subsection{" << p.begin()->_name << "}" << endl;
            out << "\\label{" << label << "}" << endl;
            out << help << endl ;
            
            out << "\\rowcolors{1}{invisiblegray}{white}" << endl;
            out << "{ \\small" << endl;
            out << "\\begin{longtable}{m{3cm}|m{11cm}}" << endl;
            
            p.PrintNodeTEX(out, "", p, -10);
            
            out << "\\end{longtable}" << endl;
            out << "}" << endl;
            break;
        }

        out << endl;
    }

    return out;
};

//{
//      Property::PrintNode(out, "", p);
      //Property::PrintNodeXML(out, "", p, "");
//      return out;
//}

static void start_hndl(void *data, const char *el, const char **attr)
{
    stack<Property *> *property_stack =
        (stack<Property *> *)XML_GetUserData((XML_Parser*)data);

    Property *cur = property_stack->top();
    Property &np = cur->add(el, "");
    
    for (int i = 0; attr[i]; i += 2)
        np.setAttribute(attr[i], attr[i + 1]);    
    
    property_stack->push(&np);
}

static void end_hndl(void *data, const char *el)
{
    stack<Property *> *property_stack =
        (stack<Property *> *)XML_GetUserData((XML_Parser*)data);
    property_stack->pop();
}

void char_hndl(void *data, const char *txt, int txtlen)
{
    stack<Property *> *property_stack =
        (stack<Property *> *)XML_GetUserData((XML_Parser*)data);

    Property *cur = property_stack->top();
    cur->value().append(txt, txtlen);
}

bool load_property_from_xml(Property &p, string filename)
{
  XML_Parser parser = XML_ParserCreate(NULL);
  if (! parser)
    throw std::runtime_error("Couldn't allocate memory for xml parser");

  XML_UseParserAsHandlerArg(parser);
  XML_SetElementHandler(parser, start_hndl, end_hndl);
  XML_SetCharacterDataHandler(parser, char_hndl);

  ifstream fl;
  fl.open(filename.c_str());
  if(!fl.is_open())
    throw std::ios_base::failure("Error on open xml file: " + filename);

  stack<Property *> pstack;
  pstack.push(&p);

  XML_SetUserData(parser, (void*)&pstack);
  while(!fl.eof()) {
    string line;
    getline(fl, line);
    line=line + "\n";
    if (! XML_Parse(parser, line.c_str(), line.length(), fl.eof()))
      throw  std::ios_base::failure(filename + ": Parse error at line " +
          boost::lexical_cast<string>(XML_GetCurrentLineNumber(parser)) + "\n" +
          XML_ErrorString(XML_GetErrorCode(parser)));
  }
  fl.close();
  return true;
}

}}
