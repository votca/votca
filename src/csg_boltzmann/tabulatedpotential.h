/* 
 * File:   tabulatedpotential.h
 * Author: ruehle
 *
 * Created on November 4, 2009, 4:53 PM
 */

#ifndef _TABULATEDPOTENTIAL_H
#define	_TABULATEDPOTENTIAL_H

#include <vector>
#include "bondedstatistics.h"
#include <votca/tools/histogram.h>

using namespace std;

class TabulatedPotential
    : public AnalysisTool
{
    public:
        TabulatedPotential();
        ~TabulatedPotential() {};

        void Register(map<string, AnalysisTool *> &lib);

        void Command(BondedStatistics &bs, string cmd, vector<string> &args);

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

