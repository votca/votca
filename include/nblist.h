/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

