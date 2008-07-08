// 
// File:   correlate.cc
// Author: ruehle
//
// Created on July 30, 2007, 10:58 AM
//

#include "correlate.h"
#include <math.h>

/**
    \todo clean implementation!!!
*/
void Correlate::CalcCorrelations(DataCollection<double>::selection *data)
{    
    size_t N;
    double xm(0), xsq(0);    
    
    N = (*data)[0].size();
     for(int i=0; i<N; i++) {
        xm += (*data)[0][i];
        xsq += (*data)[0][i]*(*data)[0][i];
    }
    xm/=(double)N;
    
    for(int v=1; v<data->size(); v++) {
        pair<string, double> p("do_names", 0);        
        double ym(0), ysq(0);
        
        for(int i=0; i<N; i++) {
            ym+=(*data)[v][i];
            ysq+=(*data)[v][i]*(*data)[v][i];
            p.second+= (*data)[v][i]*(*data)[0][i];
        }
        ym/=(double)N;
        double norm = (xsq - ((double)N)*xm*xm)*(ysq - ((double)N)*ym*ym);
        p.second = (p.second - ((double)N)*xm*ym) / sqrt(norm);
        _corr.push_back(p);
    }
}
