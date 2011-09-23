/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _QMNBList_H
#define	_QMNBList_H

#include "qmpair.h"
#include <votca/csg/pairlist.h>
#include <votca/csg/beadlist.h>
#include <votca/moo/crgunit.h>

namespace votca { namespace ctp {


using namespace votca::csg;

class QMTopology;
using namespace votca::tools;
using namespace votca::csg;

/**
 * \brief Neighbour search for crg units
 *
 * This class wraps the NBList from csg to work on CrgUnits,
 * this all looks a bit cumbersome now, but will make things
 * nice if one want's to switch in between nbsearch algorithms
 *
 * */

class QMNBList
    : public PairList<QMCrgUnit *, QMPair>
{
public:
     QMNBList(): _cutoff(0.) {};
     ~QMNBList(){
         PairList<QMCrgUnit *, QMPair>::Cleanup();
     }
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true) { Generate(list, list, do_exclusions); }

    void setCutoff(double cutoff) { _cutoff = cutoff; }
    double getCutoff() {return _cutoff; }

protected:

    bool Match(Bead *b1, Bead *b2, const vec &r, const double notused);
    double _cutoff;
    QMTopology *_father;

};

}}

#endif	/* _QMNBList_H */

