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

#ifndef _HISTOGRAMNEW_H
#define	_HISTOGRAMNEW_H

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

#include "table.h"

namespace votca { namespace tools {

using namespace std;

/**
    \brief class to generate histograms

    This class produces a histogram out of a vector of values
*/
class HistogramNew
{
    public:                
       
        /// constructor
        HistogramNew();
        /// destructor
        ~HistogramNew() {};
        
        /**
         * \brief Initialize the Histogram
         * @param min lower bound of interval
         * @param max upper bound of interval
         * @param nbins number of bins
         */void Initialize(double min, double max, int nbins);
        
        /**
          * \brief process a data point
          * \param v value of this point
          * \scale scale weighting of this point, bin of v is increased by scale instead of 1
         */
        void Process(const double &v, double scale = 1.0);
        
        /**
            \brief process a range of data using iterator interface
         */
        template<typename iterator_type>        
        void ProcessRange(const iterator_type &begin, const iterator_type &end);
    
    
        /**
         * \brief get the lower bound of the histogram intervaö
         * \return lower limit of interval
         */
        double getMin() const {return _min; }
        /**
         * \brief get the upper bound of the histogram intervaö
         * \return upper limit of interval
         */
        double getMax() const {return _max; }
        /**
         * \brief Get number of grid points
         * \return number of grid poitns
         */
        double getNBins() const {return  _nbins; }

        /**
         * \brief get the grid of histogram
         * \return step per bin
         */
        double getStep() const { return _step; }
        
        /**
         * \brief normalize the histogram that the integral is 1
         */
        void Normalize();
        
        /**
         * \brief clear all data
         */
        void Clear();
        
    
        /**
         * \brief get access to content of histogram
         * \return table object with bins in x and values in y
         */
        Table &data() { return _data; }
        
    private:        
        double _weight;
        double _min, _max;
        double _step;
        bool _periodic; 
        
        int _nbins;
        
        Table _data;
};

inline ostream& operator<<(ostream& out, HistogramNew &h)
{
    out << h.data();
    //for(int i=0; i<h.getNBins(); i++) {
    //    out << h.getMin() + h.getStep()*((double)i+0.0) << " " << h._get()[i] << endl;
    //}
    return out;
}

template<typename iterator_type>
inline void HistogramNew::ProcessRange(const iterator_type &begin, const iterator_type &end)
{
    for(iterator_type iter = begin; iter!=end; ++iter)
        Process(*iter);
}

}}

#endif	/* _HISTOGRAMNEW_H */

