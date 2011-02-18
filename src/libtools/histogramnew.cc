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

#include "histogramnew.h"

namespace votca { namespace tools {

HistogramNew::HistogramNew()
{
    _min=_max=_step=0;
    _weight = 1.;
    _periodic=false;
}

/// \todo implement this correctly
HistogramNew::HistogramNew(const HistogramNew &hist)
{
    HistogramNew();
}


void HistogramNew::Initialize(double min, double max, int nbins)
{
    _min = min; _max = max;
    _step = (_max - _min)/nbins;
    _weight = 1.;
    _data.resize(nbins);  
    _nbins = nbins;
    
    for(double v=_min, i=0; i<nbins; v+=_step,++i)
        _data.x(i)=v;
    
    _data.y()=ub::zero_vector<double>(_nbins);
    _data.yerr()=ub::zero_vector<double>(_nbins);
    _data.flags()=ub::scalar_vector<char>(_nbins, 'i');    
}

void HistogramNew::Process(const double &v, double scale)
{
    int i = (int) ((v - _min) / _step + 0.5);
    
    if (i < 0 || i >= _nbins) {
        if(!_periodic) return;
        if(i<0) i = _nbins - ((-i) % _nbins);
        else i = i % _nbins;        
    }
    _data.y(i) += _weight * scale;
} 

void HistogramNew::Normalize()
{
    double area = 0;
    
    
    area=ub::norm_1(_data.x()) * _step;
    
    _weight /= area;
    double scale = 1./area;
    
    _data.y() *= scale;    
}

void HistogramNew::Clear()
{
    _weight = 1.;
    _data.y() = ub::zero_vector<double>(_nbins);
    _data.yerr() = ub::zero_vector<double>(_nbins);
}

}}
