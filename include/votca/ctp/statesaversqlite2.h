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

#ifndef __VOTCA_MD2QM_STATESAVERSQLITE2_H
#define	__VOTCA_MD2QM_STATESAVERSQLITE2_H

#include <stdio.h>
#include <map>
// #include "qmbead.h"
//#include "qmtopology.h"
// #include "qmpair.h"
// #include "qmnblist.h"
#include "qmdatabase2.h"
#include "topology.h"

namespace votca { namespace ctp {

using namespace votca::tools;

class StateSaverSQLite2
{
public:
    StateSaverSQLite2() { };
   ~StateSaverSQLite2() { _db.Close(); }

    void Open(Topology &qmtop, const string &file);
    void Close() { _db.Close(); }
    bool NextFrame();

    void WriteFrame();
    void WriteMeta(bool update);
    void WriteMolecules(bool update);
    void WriteSegments(bool update);
    void WriteFragments(bool update);
    void WriteAtoms(bool update);
    void WritePairs(bool update);

    void ReadFrame();
    void ReadMeta(int topId);
    void ReadMolecules(int topId);
    void ReadSegments(int topId);
    void ReadFragments(int topId);
    void ReadAtoms(int topId);
    void ReadPairs(int topId);

    int  FramesInDatabase();
    Topology *getTopology() { return _qmtop; }
    bool HasTopology(Topology *top);
    
private:
    Topology       *_qmtop;
    QMDatabase2     _db;

    int             _frame;
    int             _current_frame;
    vector<int>     _frames;
    vector<int>     _topIds;

    string          _sqlfile;
    bool            _was_read;
};

}}

#endif	/* __VOTCA_MD2QM_STATESAVERSQLITE2_H */

