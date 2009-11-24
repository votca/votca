/* 
 * File:   nblistgrid.h
 * Author: ruehle
 *
 * Created on October 26, 2009, 4:21 PM
 */

#ifndef _NBLISTGRID_H
#define	_NBLISTGRID_H

#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>
#include "nblist.h"
#include <vector>

class NBListGrid
    : public NBList
{
public:
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true);

protected:
    struct cell_t {
        BeadList _beads;
        std::vector<cell_t*> _neighbours;
    };

    vec _box_a, _box_b, _box_c;
    vec _norm_a, _norm_b, _norm_c;
    int _box_Na, _box_Nb, _box_Nc;

    std::vector<cell_t> _grid;
    Topology *_top;

    void InitializeGrid(const matrix &box);
    
    cell_t &getCell(const vec &r);
    cell_t &getCell(const int &a, const int &b, const int &c);

    void TestBead(cell_t &cell, Bead *bead);
    void TestCell(cell_t &cell, Bead *bead);
};

inline NBListGrid::cell_t &NBListGrid::getCell(const int &a, const int &b, const int &c)
{
    return _grid[a + _box_Na*b + _box_Na*_box_Nb*c];
}

#endif	/* _NBLISTGRID_H */

