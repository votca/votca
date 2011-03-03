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

#include <votca/tools/tokenizer.h>
#include "readxml.h"

void ReadXML::Initialize(QMTopology* top, Property* options)
{
     _filename = options->get("options.readxml.file").as<string>();
}

bool ReadXML::EvaluateFrame(QMTopology *top)
{    
    _top = top;
    _parser.NextHandler(this, &ReadXML::ParseRoot);
    _parser.Open(_filename);
}

void ReadXML::EndEvaluate(QMTopology *top)
{
}


void ReadXML::ParseRoot(const string &el, map<string, string> &attr)
{
    if(el != "qmtop")
        throw std::runtime_error("Wrong root node in xml. Must be qmtop.");
    _parser.NextHandler(this, &ReadXML::ParseBody);
}

void ReadXML::ParseBody(const string &el, map<string, string> &attr)
{
    if(el == "frame") {
        _parser.NextHandler(this, &ReadXML::ParseFrame);
    }

    else
        throw std::runtime_error("error, unknown node: " + el);
}

void ReadXML::ParseFrame(const string &el, map<string, string> &attr)
{
    if(el == "clear_nblist") {
        _top->nblist().Cleanup();        
    }
    else if(el == "pair") {

        int first = lexical_cast<int>(attr["first"]) - 1;
        int second = lexical_cast<int>(attr["second"]) - 1;

        QMCrgUnit *crg1 = _top->GetCrgUnit(first);
        QMCrgUnit *crg2 = _top->GetCrgUnit(second);

        if(crg1->getId() > crg2->getId())
            swap(crg1, crg2);

        if(_top->nblist().FindPair(crg1, crg2))
            throw std::runtime_error("multiple definitions of pair (" 
                    + lexical_cast<string>(first+1) + ", "
                    + lexical_cast<string>(second+1) + ")");
        QMPair *pair = new QMPair(crg1, crg2, _top);

        map<string,string>::iterator iter;
        for(iter=attr.begin(); iter!=attr.end(); ++iter) {
            if(iter->first == "J") {
                Tokenizer tok(iter->second, " ");
                tok.ConvertToVector<double>(pair->Js());

            }
            else if(iter->first == "rate12") {
                pair->setRate12(lexical_cast<double>(iter->second));
            }
            else if(iter->first == "rate21") {
                pair->setRate21(lexical_cast<double>(iter->second));
            }
            else if(iter->first == "first" || iter->first == "second") {
            }
            else if(iter->first == "lambda_out") {
                pair->setLambdaOuter(lexical_cast<double>(iter->second));
            }
            else
                throw std::runtime_error("undefined property in pair: \"" + iter->first + "\"");
        }
        _top->nblist().AddPair(pair);
        
    }
    else if (el == "site" ) {
    	int site_number = lexical_cast<int>(attr["number"]) - 1;
	CrgUnit *crg = _top->GetCrgUnit(site_number);

        map<string,string>::iterator iter;
	for(iter=attr.begin(); iter!=attr.end(); ++iter) {
	    if(iter->first == "energy") {
		crg->setEnergy(lexical_cast<double>(iter->second));
	    }
	}
    }
    else throw std::runtime_error("unknown node: " + el);
    _parser.IgnoreChilds();
}
