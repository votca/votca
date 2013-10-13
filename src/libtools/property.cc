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
#include <stdexcept>
#include <stdio.h>
#include <expat.h>
#include <string.h>
#include <fstream>
#include <string>
#include <stack>
#include <iomanip>

#include <votca/tools/property.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/propertyformat.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace votca { namespace tools {

// ostream modifier defines the output format, level, indentation
const int Property::_format = std::ios_base::xalloc();  
const int Property::_output_level = std::ios_base::xalloc();  
const int Property::_output_indent = std::ios_base::xalloc();  
   
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
