/* 
 * File:   imcoefficients.cpp
 * Author: victorr
 *
 * Created on June 10, 2008, 12:05 PM
 */

#include "imccoefficients.h"
#include "topology.h"

void IMCCoefficients::clear()
{
    for(vector<double>::iterator i = _cross.begin(); i!=_cross.end();++i)
        *i=0;
    for(vector<double>::iterator i = _dist.begin(); i!=_dist.end();++i)
        *i=0;
    _c = 0;
}

void IMCCoefficients::Process(vector<int> dist)
{
    for(int i=0; i<_N; i++) {
        _dist[i]+=dist[i];
        for(int j=0; j<=i; j++)
            _cross[i*(i+1)/2 + j] += dist[i]*dist[j];
    }
    _c++;
}

void IMCCoefficients::Average()
{
    for(vector<double>::iterator i = _cross.begin(); i!=_cross.end();++i)
        *i/=(double)_c;
    for(vector<double>::iterator i = _dist.begin(); i!=_dist.end();++i)
        *i/=(double)_c;    
}

void IMCCoefficients::OutputCross(ostream &out)
{
    out << _N << endl;
    for(int i=0; i<_N; ++i) {
        for(int j=0; j<=i; ++j)
            out << _cross[i*(i+1)/2 + j] << " ";
        out << endl;
    }
}

void IMCCoefficients::OutputDist(ostream &out)
{
    out << _N << endl;
    for(vector<double>::iterator i = _dist.begin(); i!=_dist.end();++i)
        out << *i << endl;
}

