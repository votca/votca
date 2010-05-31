/*
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <stdio.h>
#include <stdexcept>
#include <expat.h>
#include "parsexml.h"

static void start_hndl(void *data, const char *el, const char **attr) {
    ParseXML *reader =
            (ParseXML*) XML_GetUserData((XML_Parser*) data);

    map<string, string> mattr;

    for (int i = 0; attr[i]; i += 2)
        mattr[attr[i]] = attr[i + 1];
    string sel = el;
    reader->StartElemHndl(sel, mattr);
}

static void end_hndl(void *data, const char *el) {
    ParseXML *reader =
            (ParseXML*) XML_GetUserData((XML_Parser*) data);
    reader->EndElemHndl(el);
}

void ParseXML::ParseIgnore(const string &el, map<string, string> &attr) {
    SetHandler(&ParseXML::ParseIgnore);
}

void ParseXML::StartElemHndl(const string &el, map<string, string> &attr) {
    (this->*_handler)(el, attr);
}

void ParseXML::EndElemHndl(const string &el) {
    _stack_handler.pop();
    _handler = _stack_handler.top();
}

void ParseXML::SetHandler(ElemHandler_t handler) {
    _handler = handler;
    _stack_handler.push(handler);
}

