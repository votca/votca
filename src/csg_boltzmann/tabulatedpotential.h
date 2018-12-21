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

#include <vector>
#include "bondedstatistics.h"
#include "analysistool.h"
//#include <votca/tools/histogramnew.h>
#include <votca/tools/histogram.h>

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

        /**
         * \brief Returns the temperature used during the bolzmann inversion 
         *
         * \return temperature in kelvin
         **/
        double getTemperature() const;

        /**
         * \brief Method returns the number of smoothing iterations used on the
         * data
         *
         * The first integer is the number of smoothing iterations before 
         * boltzmann inversion is done, and the second interger is the number
         * of iterations done on the boltzmann inverted histogram. 
         *
         * \return pair of integers containing the number of smoothing 
         * iteratiosn before and after boltzmann inversion.
         **/
        std::pair<int,int> getSmoothIterations() const;
    private:
        bool SetOption_(Histogram::options_t &op, const vector<string> &args);
        bool SetOption_(const vector<string> &args);

        /**
         * \brief Smooths a vector of doubles
         *
         * This function uses a weighted smoothing algorithm using 5 data 
         * points, the weights are applied as 1:2:3:2:1 where the middle point
         * becomes the average of these weighted values.
         *
         * \param[in,out] data vector of doubles that is smoothed
         * \param[in] periodic boolean determining if the smoothing will use
         * periodic boundary conditions
         **/
        void Smooth_(vector<double> &data, bool bPeriodic);
        void BoltzmannInvert_(vector<double> &data);
        void CalcForce_(vector<double> &u, vector<double> &du, double dx, bool bPeriodic);

        Histogram::options_t _tab_options;
        Histogram::options_t _hist_options;

        /// How many times the data is smoothed before the histogram is 
        /// boltzmann inverted.
        int _tab_smooth1;
        /// How many times the data is smoothed after the histogram is boltzmann
        /// inverted.
        int _tab_smooth2;
        /// Temperature in units of Kelvin
        double _Temperature;

//        HistogramNew _hist;
};

#endif	/* _TABULATEDPOTENTIAL_H */

