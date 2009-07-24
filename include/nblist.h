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
    NBList() : _excl(0) {}
    

    void Generate(BeadList &list1, BeadList &list2, ExclusionList *ExcList=0);
    void Generate(BeadList &list, ExclusionList *ExcList) { Generate(list, list, ExcList); }
    
    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() { return _cutoff; }
   
private:
    double _cutoff;
    ExclusionList *_excl;
    
    bool Match(Bead *bead1, Bead *bead2, const vec &r, ExclusionList *ExcList=0);
};

#endif	/* _NBLIST_H */

