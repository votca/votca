/// \addtogroup csg
///@{
/* 
 * File:   beadlist.h
 * Author: ruehle
 *
 * Created on July 17, 2008, 5:14 PM
 */

#ifndef _BEADLIST_H
#define	_BEADLIST_H

#include <string>

using namespace std;

class BeadList
{
public:
    BeadList() {}
    ~BeadList() {}
    
    int CreateList(string select);
    
private:
};

#endif	/* _BEADLIST_H */

/// @}
