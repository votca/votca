/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _correlate_H
#define	_correlate_H


#include <vector>
#include <iostream>
#include "datacollection.h"

namespace votca { namespace tools {

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

}}

#endif	/* _correlate_H */

