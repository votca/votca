/* 
 * File:   histogramnew.cc
 * Author: ruehle
 *
 * Created on March 11, 2009, 4:29 PM
 */

#include "histogramnew.h"

HistogramNew::HistogramNew()
{
    _min=_max=_step=0;
    _weight = 1.;
    _periodic=false;
}

void HistogramNew::Initialize(double min, double max, int nbins)
{
    _min = min; _max = max;
    _step = (_max - _min)/nbins;
    _weight = 1.;
    _data.resize(nbins);  
    _nbins = nbins;
    
    int i;
    for(double v=_min, i=0; i<nbins; v+=_step,++i)
        _data.x(i)=v;
    
    _data.y()=ub::zero_vector<double>(_nbins);    
    _data.flags()=ub::scalar_vector<char>(_nbins, 'i');    
}

void HistogramNew::Process(double &v)
{
    int i = (int) ((v - _min) / _step + 0.5);
    
    if (i < 0 || i >= _nbins) {
        if(!_periodic) return;
        if(i<0) i = _nbins - ((-i) % _nbins);
        else i = i % _nbins;        
    }
    _data.y(i) += _weight;
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
}