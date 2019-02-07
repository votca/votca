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

#ifndef _TABULATEDPOTENTIAL_H
#define	_TABULATEDPOTENTIAL_H

#include "analysistool.h"
#include "bondedstatistics.h"
#include <votca/tools/histogram.h>
#include <vector>

using namespace std;
using namespace votca::tools;
using namespace votca::csg;

class TabulatedPotential
    : public AnalysisTool
{
    public:
        TabulatedPotential();
        ~TabulatedPotential() {};

        void Register(map<string, AnalysisTool *> &lib);

        void Command(BondedStatistics &bs, string cmd, vector<string> &args);
        void Help(string cmd, vector<string> &args);

        void WriteHistogram(BondedStatistics &bs, vector<string> &args);
        void WritePotential(BondedStatistics &bs, vector<string> &args);

    private:
        bool SetOption(Histogram::options_t &op, const vector<string> &args);
        void Smooth(vector<double> &data, bool bPeriodic);
        void BoltzmannInvert(vector<double> &data, double T);
        void CalcForce(vector<double> &u, vector<double> &du, double dx, bool bPeriodic);

        Histogram::options_t _tab_options;
        Histogram::options_t _hist_options;
        int _tab_smooth1;
        int _tab_smooth2;
        double _T;
};

#endif	/* _TABULATEDPOTENTIAL_H */

