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
    NBList();
    
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true) { Generate(list, list, do_exclusions); }
    
    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() { return _cutoff; }

    /// functon to use a user defined pair type
    template<typename pair_type>
    void setPairType();

    
    /// typedef for a user match function, return true if bead should be added
    typedef bool (*match_function_t)(Bead *, Bead *, const vec &r);

    /// set user match function
    void setMatchFunction(match_function_t match_function);

    /// match function that always matches
    static bool match_always(Bead *b1, Bead *b2, const vec &r) { return true; }

protected:
    double _cutoff;
    bool _do_exclusions;

    /// policy function to create new bead types
    template<typename pair_type>
    static BeadPair *beadpair_create_policy(Bead *bead1, Bead *bead2, const vec &r)
    {
        return new pair_type(bead1, bead2, r);
    }

    typedef BeadPair* (*pair_creator_t)(Bead *bead1, Bead *bead2, const vec &r);
    /// the current bead pair creator function
    pair_creator_t _pair_creator;

    //    typedef T* (*creator_t)();

    match_function_t _match_function;
};

template<typename pair_type>
void NBList::setPairType()
{
    _pair_creator = NBList::beadpair_create_policy<pair_type>;
}

inline void NBList::setMatchFunction(match_function_t match_function)
{
    _match_function = match_function;
}

#endif	/* _NBLIST_H */

