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

/**
 * \brief Tabulated Potential calculates histograms of bead interactions
 *
 * There are two histograms that can be created:
 *
 * * the first is a potential distribution table specified with the 'tab' 
 * keyword (for table). This uses boltzmann inversion to calculate the 
 * potential from a selection of bead interactions.
 *
 * * the second is a histogram of the beads interactions specified with the 
 * 'hist' keyword. 
 *
 * As an example let us assume we already have a BondedStatistics object named
 * 'bonded_statistics' with the correct data. We can write a histogram of the 
 * interactions using the following commands:
 *
 *     TabulatedPotential tabulated_potential;
 *
 *     map<string,AnalysisTool *> commands;
 *     tabulated_potential.Register(commands);
 *
 *     vector<string> arguments{"set","n","11"};
 *
 *     tabulated_potential.Command(bonded_statistics,"hist",arguments);
 *
 *     vector<string> interactions{"file1.txt","bond1","bond2","bond3",...};
 *     tabulated_potential.Command(bonded_statistics,"hist",interactions);
 *
 * The above example creates the object it then proceeds to set one of the
 * properties of the histogram. In this case the number of data points is
 * specified with the keyword 'n', the 11 indicates we want the histogram to
 * consist of 11 bins. To actually create the histogram the Command method is 
 * called a second time but because the 'set' keyword is not the first element
 * in the vector of strings the tabulated potential object assumes it is a file
 * name were the histogram will be written. The following labels ("bond1",
 * "bond2", etc...) specify the interactions to be sorted into the histogram. 
 * Running this code should produce a file named file1.txt containing two 
 * columns. 
 *
 * column1 = the bin edge, column2 = the contents of the histogram
 *
 * If the table potential is printed instead of the historgram of interactions
 * than 3 columns are printed e.g. in the example below:
 *
 *     TabulatedPotential tabulated_potential;
 *
 *     map<string,AnalysisTool *> commands;
 *     tabulated_potential.Register(commands);
 *
 *     vector<string> arguments{"set","T","275"};
 *     vector<string> arguments2{"set","smooth_pot","1"};
 *     vector<string> arguments3{"set","periodic","1"};
 *
 *     tabulated_potential.Command(bonded_statistics,"tab",arguments);
 *     tabulated_potential.Command(bonded_statistics,"tab",arguments1);
 *     tabulated_potential.Command(bonded_statistics,"tab",arguments2);
 *
 *     vector<string> interactions{"file2.txt","bond1","bond2","bond3",...};
 *     tabulated_potential.Command(bonded_statistics,"tab",interactions);
 *
 * Here, notice we are using the 'tab' keyword to indicate it is for the 
 * tabulated potential. We are also setting several properities. The temperature
 * is changed from the default 300 Kelvin to 275 Kelvin. Smooth potential is set
 * to 1 this has the affect of smoothing the potential after boltzmann inversion
 * is done. If a value greater than 1 were used it would smooth the data several
 * times upto the number specified. And the final setting was to make the 
 * tabulated potential periodic. The value "1" is converted to a boolean which
 * is interpreted as true. 
 *
 * Finally, the data is printed to a file called "file2.txt" which contains 
 * three columns:
 * 
 * column1 = the bin edge, column2 = potential, column3 = the force
 **/
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

