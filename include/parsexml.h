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

#ifndef _PARSEXML_H
#define _PARSEXML_H

#include <string>
#include <stack>
#include <map>
#include <list>

using namespace std;

class ParseXML {
public:
    /// start element callback for xml parser
    void StartElemHndl(const string &el, map<string, string> &attr);
    /// end element callback for xml parser
    void EndElemHndl(const string &el);
    // virtual void ParseRoot(const string &el, map<string, string> &attr);
    void ParseIgnore(const string &el, map<string, string> &attr);

    typedef void (ParseXML::*ElemHandler_t)(const string &, map<string, string> &);
    // void SetHandler(ElemHandler_t handler);
    void ParseFrame(const string &el, map<string, string> &attr);

private:
    void SetHandler(ElemHandler_t handler);

    stack<ElemHandler_t> _stack_handler;
    ElemHandler_t _handler;
};

#endif