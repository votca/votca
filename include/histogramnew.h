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
        
        void Initialize(double min, double max, int nbins);        
        
        /**
            process data and generate histogram
         */
        void Process(double &v, double scale = 1.0);
        
        /**
            process a range of data
         */
        template<typename iterator_type>        
        void ProcessRange(const iterator_type &begin, const iterator_type &end);
    
    
        /// returns the minimum value
        double getMin() const {return _min; }
        /// return the maximum value
        double getMax() const {return _max; }
        /// return the number of grid points
        double getNBins() const {return  _nbins; }

        double getStep() const { return _step; }             
        
        /// normalize the histogram
        void Normalize();
        
        void Clear();
        
    
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

