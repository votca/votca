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
 * File:   stdanalysis.h
 * Author: ruehle
 *
 * Created on November 4, 2009, 4:53 PM
 */

#ifndef _STDANALYSIS_H
#define	_STDANALYSIS_H

#include "bondedstatistics.h"
#include <map>

using namespace std;

class StdAnalysis
    : public AnalysisTool
{
    public:
        StdAnalysis() {};
        ~StdAnalysis() {};

        void Register(map<string, AnalysisTool *> &lib);

        void Command(BondedStatistics &bs, string cmd, vector<string> &args);
        void Help(string cmd, vector<string> &args);

        void WriteValues(BondedStatistics &bs, vector<string> &args);
        void WriteCorrelations(BondedStatistics &bs, vector<string> &args);
        void WriteAutocorrelation(BondedStatistics &bs, vector<string> &args);
    private:
};


#endif	/* _STDANALYSIS_H */

