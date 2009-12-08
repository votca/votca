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
/* 
 * File:   analysistool.h
 * Author: ruehle
 *
 * Created on August 2, 2007, 3:12 PM
 */

#ifndef _analasystool_H
#define	_analasystool_H

#include <map>
#include <string>
#include "cgengine.h"
#include "bondedstatistics.h"

using namespace std;

/**
    \brief base class for all analasys tools

    This is the base class for all analasys tool. 
    \todo do option functions!!!
*/
class AnalysisTool
{
public:
    AnalysisTool() {}
    virtual ~AnalysisTool() {}
    
    virtual void Register(map<string, AnalysisTool *> &lib) {}
    virtual void Command(BondedStatistics &bs, string cmd, vector<string> &args) {};
    virtual void Help(string cmd, vector<string> &args) {};
    
private:
//    map<string, string> _options;
};

#endif	/* _analasystool_H */

