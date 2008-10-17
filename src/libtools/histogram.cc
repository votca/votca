// 
// File:   histogram.cc
// Author: ruehle
//
// Created on April 26, 2007, 4:53 PM
//

#include <limits>
#include <math.h>
#include "histogram.h"

Histogram::Histogram()
    : _min(0), _max(0)
{}

Histogram::Histogram(options_t op)
    : _options(op)
{    
}

Histogram::~Histogram()
{
}

void Histogram::ProcessData(DataCollection<double>::selection *data)
{
    DataCollection<double>::selection::iterator array;
    DataCollection<double>::array::iterator iter;
    int ii;
    int ndata = 0;
    
    _pdf.assign(_options._n, 0);
        
    if(_options._auto_interval) {
        _min = numeric_limits<double>::max();
        _max = numeric_limits<double>::min();
        _options._extend_interval = true;
    }
    else {
        _min = _options._min;
        _max = _options._max;
    }
    
    for(array = data->begin(); array!=data->end(); ++array) {
        ndata+=(*array)->size();
        if(_options._extend_interval || _options._auto_interval) {
            for(iter=(*array)->begin(); iter!=(*array)->end(); ++iter) {
                _min = min(*iter, _min);
                _max = max(*iter, _max);            
            }
        }
    }
    
    // make that the highes value fits into interval
    //if(_options._auto_interval || _max!=_options._max)
    //    _max = _max + 0.5*(_max - _min)/(double)(_options._n);
    
    _interval = (_max - _min)/(double)(_options._n-1);

    double v = 1./(double)ndata/_interval;
    for(array = data->begin(); array!=data->end(); ++array) {
        for(iter=(*array)->begin(); iter!=(*array)->end(); ++iter) {
            ii = (int)( (*iter - _min) / _interval + 0.5); // the interval should be centered around the sampling point              
            if(ii< 0 || ii >= _options._n) {
                if(_options._periodic) {
                    while(ii<0) ii+=_options._n;
                    ii = ii % _options._n;
                }
                else { continue; } //cout << "[histogram.cc]: out of bounds" << endl; continue;}
            }
            _pdf[ii]+= v;
        }
    }        
    
    //cout << _pdf.size() << " " << _options._periodic << endl;
    if(_options._scale == "bond") {
        for(int i=0; i<_pdf.size(); ++i) {
            double r = _min + _interval*(double)i;
            if(abs(r) < 1e-10) {
                r = _min + _interval*(double)(i+1);
                _pdf[i] = _pdf[i+1];
            }
            _pdf[i] /= (r*r);
        }
    }
    else if(_options._scale == "angle") {
        for(int i=0; i<_pdf.size(); ++i) {
            double alpha = _min + _interval*(double)i;
            double sa = sin(alpha); 
            if(abs(sa) < 1e-5) {
                if(i<_pdf.size()-1) {
                    alpha = _min + _interval*(double)(i+1);
                    _pdf[i]=_pdf[i+1]/sin(alpha);
                }
                else {
                    _pdf[i]=_pdf[i-1];
                }
                
            }
            else _pdf[i] /= sa;
        }        
    }
    
    if(_options._periodic) {
        _pdf[0] = (_pdf[0] + _pdf[_options._n-1]);
        _pdf[_options._n-1] = _pdf[0];
    }
}

