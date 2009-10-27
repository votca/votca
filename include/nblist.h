/* 
 * File:   nblist.h
 * Author: ruehle
 *
 * Created on March 5, 2009, 3:07 PM
 */

#ifndef _NBLIST_H
#define	_NBLIST_H

#include "beadlist.h"
#include "beadpair.h"
#include "pairlist.h"
#include "exclusionlist.h"

/**
 * \brief new implementation of neighbourlist, will substitute Neighbourlist
 * 
 */
class NBList 
    : public PairList<Bead*, BeadPair>
{
public:
    NBList() : _do_exclusions(false) {}
    
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true) { Generate(list, list, do_exclusions); }
    
    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() { return _cutoff; }
   
protected:
    double _cutoff;
    bool _do_exclusions;
    bool Match(Bead *bead1, Bead *bead2, const vec &r);
};

#endif	/* _NBLIST_H */

