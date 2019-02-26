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

#include <boost/lexical_cast.hpp>
#include <expat.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <votca/tools/parsexml.h>

namespace votca {
namespace tools {

void start_hndl(void *data, const char *el, const char **attr) {
  ParseXML *reader = (ParseXML *)XML_GetUserData((XML_Parser *)data);

  map<string, string> mattr;

  for (int i = 0; attr[i]; i += 2) mattr[attr[i]] = attr[i + 1];
  string sel = el;
  reader->StartElemHndl(sel, mattr);
}

void end_hndl(void *data, const char *el) {
  ParseXML *reader = (ParseXML *)XML_GetUserData((XML_Parser *)data);
  reader->EndElemHndl(el);
}

void ParseXML::Open(const string &filename) {
  XML_Parser parser = XML_ParserCreate(NULL);
  if (!parser)
    throw std::runtime_error("Couldn't allocate memory for xml parser");

  XML_UseParserAsHandlerArg(parser);
  XML_SetElementHandler(parser, start_hndl, end_hndl);
  //    XML_SetCharacterDataHandler(parser, char_hndl);

  ifstream fl;
  fl.open(filename.c_str());
  if (!fl.is_open())
    throw std::ios_base::failure("Error on open xml file: " + filename);

  XML_SetUserData(parser, (void *)this);
  while (!fl.eof()) {
    string line;
    getline(fl, line);
    line = line + "\n";
    if (!XML_Parse(parser, line.c_str(), line.length(), fl.eof()))
      throw std::ios_base::failure(
          filename + ": Parse error in " + filename + " at line " +
          boost::lexical_cast<string>(XML_GetCurrentLineNumber(parser)) + "\n" +
          XML_ErrorString(XML_GetErrorCode(parser)));
  }
  fl.close();
}

void ParseXML::ParseIgnore(const string &el, map<string, string> &attr) {
  NextHandler(this, &ParseXML::ParseIgnore);
}

void ParseXML::StartElemHndl(const string &el, map<string, string> &attr) {
  (*_handler)(el, attr);
}

void ParseXML::EndElemHndl(const string &el) {
  delete _handler;
  _stack_handler.pop();
  _handler = _stack_handler.top();
}

}  // namespace tools
}  // namespace votca
