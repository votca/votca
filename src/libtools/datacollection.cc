// 
// File:   datacollection.cc
// Author: ruehle
//
// Created on May 11, 2007, 3:18 PM
//

#include <sstream>
#include "datacollection.h"

ostream& operator<<(ostream& out, DataCollection<double>::selection &sel)
{
    DataCollection<double>::selection::iterator iter;
    if(sel.empty()) {
        out << "-- empty selection --" << endl;
        return out;
    }
    
    stringstream s;
    int written;
    for(size_t i=0; ; ++i) {
        s.clear();
        s.str("");
        s.setf(ios::scientific);
        written = 0;
        for(size_t j=0; j<sel.size(); j++) {      
            if(i >= sel[j].size()) {
                s << " -";
                continue;
            }
            written++;
            s << " " << (double)sel[j][i];
            
        }
        if(written == 0) return out;
        out << i << s.str() << endl;        
    }
    return out;
}
