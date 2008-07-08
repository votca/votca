/* 
 * File:   correlate.h
 * Author: ruehle
 *
 * Created on July 30, 2007, 10:54 AM
 */

#ifndef _correlate_H
#define	_correlate_H


#include <vector>
#include <iostream>
#include "datacollection.h"

using namespace std;

/**
    \brief class to calculate correlations of values

*/
class Correlate
{
    public:
        /// constructor
        Correlate() {};
        /// destructor
        ~Correlate() {};
                
        /**
            calculate the correlation of the first row in selection with all the other
               
         */
        void CalcCorrelations(DataCollection<double>::selection *data);

        vector< pair<string,double> > &getData() { return _corr; }
    private:
        vector< pair<string,double> > _corr;
};

inline ostream& operator<<(ostream& out, Correlate &c)
{
    vector< pair<string,double> > &data = c.getData();
    for(size_t i=0; i<data.size(); i++) {
        out << data[i].second << endl;
    }
    return out;
}

#endif	/* _correlate_H */

