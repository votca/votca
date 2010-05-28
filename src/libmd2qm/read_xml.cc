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

#include "read_xml.h"

bool ReadXML::EvaluateFrame(int nr, int nframes, QMTopology *top){
    _in_int.open("integrals.xml");
    QMNBList &nblist = top->nblist();
    nblist.Cleanup(); /// delete current nbl

}

void ReadXML::ParseRoot(const string &el, map<string, string> &attr)
{
    if(el == "qmtop") {
        SetHandler(&ParseXML::ParseFrame);
    }
    else {
        throw std::runtime_error("Wrong root node in xml integrals.xml. Must be qmtop.");
    }
}

void ReadXML::ParsePair(const string &el, map<string, string> &attr, QMNBList &nblist){
    if(el == "pair"){
        string first = attr["1st"];
        string second = attr["2nd"];
        string J_0 = attr["J_0"];
        vector <double> Js;
        Js.push_back(boost::lexical_cast<double>(J_0));
    }
    else{
        throw std::runtime_error("Wrong  subnode of frame in xml integrals.xml. Must be pair.");
    }
}