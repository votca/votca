/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#ifndef _QMNBList_H
#define	_QMNBList_H


#include <stdlib.h>
#include <votca/csg/pairlist.h>
#include <votca/ctp/qmpair.h>

namespace CSG = votca::csg;


namespace votca { namespace ctp {

class Topology;


class QMNBList : public CSG::PairList< Segment*, QMPair >
{
public:
    
    
    class SuperExchangeType
    {
    public:
        
        // TODO finish
        SuperExchangeType(string initString) : _initString(initString) { ; }
        
        string _initString;
        string donor;
        string acceptor;
        list<string> bridges;

        bool isOfBridge(string segment_type ) {
            std::list<string>::iterator findIter = std::find(bridges.begin(), bridges.end(), segment_type);
            return findIter != bridges.end();
        };

        bool isOfDonorAcceptor ( string segment_type ) {
            return segment_type == donor || segment_type == acceptor ;
        }
    };

    QMNBList() : _top(NULL), _cutoff(0) { };
    QMNBList(Topology* top) : _top(top), _cutoff(0) { };
   ~QMNBList() { 
       CSG::PairList<Segment*, QMPair>::Cleanup();       
       // cleanup the list of superexchange pairs
       for ( std::list<SuperExchangeType*>::iterator it = _superexchange.begin() ; it != _superexchange.end(); it++  ) {
           delete *it;
       }
   }
    
    void GenerateSuperExchange();
    
    void AddSuperExchangeType(string type) { _superexchange.push_back(new SuperExchangeType(type)); }
    void setSuperExchangeTypes(list<SuperExchangeType*> types) { _superexchange = types; }
    list<SuperExchangeType*> &getSuperExchangeTypes() { return _superexchange; }

    void    setCutoff(double cutoff) { _cutoff = cutoff; }
    double  getCutoff() { return _cutoff; }

    QMPair *Add(Segment* seg1, Segment* seg2);

    void PrintInfo(FILE *out);

protected:
    
    double      _cutoff;
    Topology   *_top;
    list<SuperExchangeType*> _superexchange;
};












}}


#endif	/* _QMNBList_H */

