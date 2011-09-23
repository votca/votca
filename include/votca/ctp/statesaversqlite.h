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

#ifndef __VOTCA_MD2QM_STATESAVERSQLITE_H
#define	__VOTCA_MD2QM_STATESAVERSQLITE_H

#include <stdio.h>
#include <map>
#include "qmbead.h"
#include "qmtopology.h"
#include "qmpair.h"
#include "qmnblist.h"
#include "qmdatabase.h"

namespace votca { namespace ctp {

using namespace votca::tools;

class StateSaverSQLite
{
public:
    StateSaverSQLite();
    ~StateSaverSQLite();

    void Open(QMTopology & qmtop,const string &file);
    void Close();
    void WriteFrame();

    bool NextFrame();

    int FramesInDatabase() { return _frames.size(); }
    
private:
    int _frame;
    QMDatabase _db;
    
    void WriteMolecules(int frameid);
    void WriteConjugatedSegments(int frameid);
    void WriteBeads(int frameid);
    void WritePairs(int frameid);
    void WriteIntegrals(QMPair *pair);

    void ReadFrame(void);
    void ReadMolecules(void);
    void ReadConjugatedSegments(void);
    void ReadBeads(void);
    void ReadPairs(void);
    void ReadIntegrals();

    template<typename T>
    void WriteCustomProperties(int object_id, std::map<string, T> &properties,
        string table, const string field_objectid, const string field_key="key", const string field_value="value");

    template<typename T>
    void ReadCustomProperties(int object_id, std::map<string, T> &properties,
        string table, const string field_objectid, const string field_key="key", const string field_value="value");

    QMTopology *_qmtop;

    vector<int> _frames;
    int _current_frame;

    map<int,int> _conjseg_id_map;
    bool _was_read;
};

}}

#endif	/* __VOTCA_MD2QM_StateSaverSQLite_H */

