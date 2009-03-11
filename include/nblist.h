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

/**
 * \brief new implementation of neighbourlist, will substitute Neighbourlist
 * 
 */

class NBList 
    : public PairList<Bead*, BeadPair>
{
public:
    NBList() {}
    
    void Generate(BeadList &list1, BeadList &list2);
    void Generate(BeadList &list) { Generate(list, list); }
    
    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() { return _cutoff; }
   
private:
    double _cutoff;
    
    bool Match(Bead *bead1, Bead *bead2, const vec &r);
};

#endif	/* _NBLIST_H */

