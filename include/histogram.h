/* 
 * File:   histogram.h
 * Author: ruehle
 *
 * Created on April 26, 2007, 4:48 PM
 */

#ifndef _histogram_H
#define	_histogram_H

#include <iostream>
#include <vector>
#include "datacollection.h"
#include <limits>
#include <cmath>

using namespace std;

/**
    \brief class to generate histograms

    This class produces a histogram out of a vector of values
*/
class Histogram
{
    public:        
        
        struct options_t;
        
        /// constructor
        Histogram(options_t op);
        Histogram();
        /// destructor
        ~Histogram();
        
        /**
            process data and generate histogram
         */
        void ProcessData(DataCollection<double>::selection *data);
        
        /// returns the minimum value
        double getMin() const {return _min; }
        /// return the maximum value
        double getMax() const {return _max; }
        /// return the number of grid points
        double getN() const {return _options._n; }
        vector<double> &getPdf() {return _pdf; }
        double getInterval() const { return _interval; }             
        
        struct options_t {            
            int _n;
            bool _auto_interval;
            bool _extend_interval;
            double _min, _max;
            bool _periodic;
            string _scale;
            
            options_t() {
                _n=101;
                _auto_interval = true;
                _extend_interval = false;
                _min = 0.; _max = 1.;
                _periodic = false;
                _scale = "no";
            }
        };          

    private:        
        vector<double> _pdf;
        double _min, _max;
        double _interval;
        
        options_t _options;
};

inline ostream& operator<<(ostream& out, Histogram &h)
{
    for(int i=0; i<h.getN(); i++) {
        out << h.getMin() + h.getInterval()*((double)i+0.0) << " " << h.getPdf()[i] << endl;
    }
    return out;
}

#endif	/* _histogram_H */

