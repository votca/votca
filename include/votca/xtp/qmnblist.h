/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
/// For earlier commit history see ctp commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_QMNBList_H
#define	VOTCA_XTP_QMNBList_H

#include <string>
#include <vector>
#include <list>

#include <votca/tools/tokenizer.h>
#include <votca/csg/pairlist.h>
#include <votca/xtp/qmpair.h>

namespace votca { namespace xtp {

class Topology;

class QMNBList : public csg::PairList< Segment*, QMPair >
{
public:
  
    QMNBList() : _top(NULL), _cutoff(0) { };
    QMNBList(Topology* top) : _top(top), _cutoff(0) { };
   ~QMNBList() { 
       csg::PairList<Segment*, QMPair>::Cleanup();       
   }
    
   
    void    setCutoff(double cutoff) { _cutoff = cutoff; }
    double  getCutoff() { return _cutoff; }

    QMPair *Add(Segment* seg1, Segment* seg2,bool safe=true);
    
    void AddQMNBlist(QMNBList &temp);
    
protected:
    
    Topology   *_top;
    double      _cutoff;
};

}}


#endif	// VOTCA_XTP_QMNBLIST_H 

