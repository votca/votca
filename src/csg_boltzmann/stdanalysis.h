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

        void WriteValues(BondedStatistics &bs, vector<string> &args);
        void WriteCorrelations(BondedStatistics &bs, vector<string> &args);
        void WriteAutocorrelation(BondedStatistics &bs, vector<string> &args);
    private:
};


#endif	/* _STDANALYSIS_H */

