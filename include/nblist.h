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

namespace votca { namespace csg {
using namespace votca::tools;

/**
 * \brief new implementation of neighbourlist, will substitute Neighbourlist
 * 
 */
class NBList 
    : public PairList<Bead*, BeadPair>
{
public:
    NBList();
    virtual ~NBList();    

    virtual void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    virtual void Generate(BeadList &list, bool do_exclusions = true) { Generate(list, list, do_exclusions); }
    
    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() { return _cutoff; }

    /// functon to use a user defined pair type
    template<typename pair_type>
    void setPairType();

    template<typename T>
    void SetMatchFunction(T *object, bool (T::*fkt)(Bead *, Bead *, const vec &));

    void SetMatchFunction(bool (*fkt)(Bead *, Bead *, const vec &));

    /// match function that always matches
    static bool match_always(Bead *b1, Bead *b2, const vec &r) { return true; }

protected:
    double _cutoff;
    bool _do_exclusions;

    /// policy function to create new bead types
    template<typename pair_type>
    static BeadPair *beadpair_create_policy(Bead *bead1, Bead *bead2, const vec &r)
    {
        return dynamic_cast<BeadPair*>(new pair_type(bead1, bead2, r));
    }

    typedef BeadPair* (*pair_creator_t)(Bead *bead1, Bead *bead2, const vec &r);
    /// the current bead pair creator function
    pair_creator_t _pair_creator;

    //    typedef T* (*creator_t)();

protected:
    // callback stuff
        class Functor {
    public:
        Functor() {}
        virtual bool operator()(Bead *, Bead *, const vec &) = 0;
    };

    template<typename T>
    class FunctorMember : public Functor {
    public:
        typedef bool (T::*fkt_t)(Bead *, Bead *, const vec &);

        FunctorMember(T* cls, fkt_t fkt) : _cls(cls), _fkt(fkt) {}

        bool operator()(Bead *b1, Bead *b2, const vec &r) {
            return (_cls->*_fkt)(b1, b2, r);
        }

    private:
        T* _cls;
        fkt_t _fkt;
    };

    class FunctorNonMember : public Functor {
    public:
        typedef bool (*fkt_t)(Bead *, Bead *, const vec &);
        FunctorNonMember(fkt_t fkt) : _fkt(fkt) {}

        bool operator()(Bead *b1, Bead *b2, const vec &r) {
            return (*_fkt)(b1, b2, r);
        }

    private:
        fkt_t _fkt;
    };

    Functor * _match_function;

};

template<typename pair_type>
void NBList::setPairType()
{
    _pair_creator = NBList::beadpair_create_policy<pair_type>;
}

template<typename T>
inline void NBList::SetMatchFunction(T *object, bool (T::*fkt)(Bead *, Bead *, const vec &))
{
    if(_match_function)
        delete _match_function;
    _match_function = dynamic_cast<Functor*>(new FunctorMember<T>(object, fkt));
}

inline void NBList::SetMatchFunction(bool (*fkt)(Bead *, Bead *, const vec &))
{
    if(_match_function)
        delete _match_function;
    _match_function = dynamic_cast<Functor*>(new FunctorNonMember(fkt));
}

}}

#endif	/* _NBLIST_H */

