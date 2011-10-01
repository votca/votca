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

#ifndef _EHISTOGRAM_H
#define	_EHISTOGRAM_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/paircalculator.h>
#include <votca/tools/histogramnew.h>

namespace votca { namespace ctp {

/**
	\brief Histogram of site energy differences of neighbor list pairs

Callname: ehistogram
*/
class Ehistogram : public PairCalculator
{
public:
    Ehistogram() {};
    ~Ehistogram() {};

//    const char *Description() { return "Histogram of site energy differences of neighbor list pairs"; }

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

private:
    HistogramNew histogram;
    double _min, _max;
    int _nbins;
    string _outfile;
};

inline void Ehistogram::Initialize(QMTopology *top, Property *options){
    _min = options->get("options.ehistogram.min").as<double>();
    _max = options->get("options.ehistogram.max").as<double>();
    _nbins = options->get("options.ehistogram.nbins").as<int>();
    _outfile = options->get("options.ehistogram.file").as<string>();

    histogram.Initialize(_min,_max,_nbins);
}

inline void Ehistogram::EndEvaluate(QMTopology *top){
    cout << "Writing Distribution of site energy differences to file " << _outfile << endl;
    histogram.data().Save(_outfile);
}

inline void Ehistogram::EvaluatePair(QMTopology *top, QMPair *pair){
    double energy1 = pair->first->getTotalEnergy();
    double energy2 = pair->second->getTotalEnergy();
    histogram.Process( energy1-energy2 );
}

}}

#endif	/* _EHISTOGRAM_H */
