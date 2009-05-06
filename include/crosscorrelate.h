// 
// File:   crosscorrelate.h
// Author: ruehle
//
// Created on May 21, 2007, 5:21 PM
//

#ifndef _crosscorrelate_H
#define	_crosscorrelate_H

#include <vector>
#include <iostream>
#include "datacollection.h"

using namespace std;

/**
    \brief class to calculate cross correlkations and autocorrelations

    This class produces a histogram out of a vector of values

    \todo implementation
*/
class CrossCorrelate
{
    public:
        /// constructor
        CrossCorrelate() {};
        /// destructor
        ~CrossCorrelate() {};
        
        /**
            calculate the cross correlation
         */
        //void CrossCorrelate(DataCollection<double>::selection *data1, 
        //    DataCollection<double>::selection *data2, bool average = false);
        
        /**
            calculate the auto correlation
               
         */
        void AutoCorrelate(DataCollection<double>::selection *data, bool average = false);

        // Calculates Fourier trafo and then auto correlation
        void AutoFourier(vector <double>& ivec);
        
        // Calculates Discrete Cosine trafo and then auto correlation
        void AutoCosine(vector <double>& ivec);
        
        // Calculates auto correlation via two Fourier trafos
        void AutoCorr(vector <double>& ivec);
        
        vector<double> &getData() { return _corrfunc; }
    private:
        vector<double> _corrfunc;
};

inline ostream& operator<<(ostream& out, CrossCorrelate &c)
{
    vector<double> &data = c.getData();
    for(size_t i=0; i<data.size(); i++) {
        out << i << " " << c.getData()[i] << endl;
    }
    return out;
}

#endif	/* _crosscorrelate_H */

