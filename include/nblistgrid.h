/* 
 * File:   nblistgrid.h
 * Author: ruehle
 *
 * Created on October 26, 2009, 4:21 PM
 */

#ifndef _NBLISTGRID_H
#define	_NBLISTGRID_H

#include "tools/matrix.h"
#include "nblist.h"

class NBListGrid
    : public NBList
{
public:
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true) { Generate(list, list, do_exclusions); }

protected:
    struct cell_t {
        BeadList _beads;
    };

    
    void InitializeGrid(const matrix &box);

};

#endif	/* _NBLISTGRID_H */

