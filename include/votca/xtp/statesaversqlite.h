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


#ifndef VOTCA_MD2QM_STATE_SAVER_SQLITE_H
#define	VOTCA_MD2QM_STATE_SAVER_SQLITE_H

#include <stdio.h>
#include <map>
#include <votca/xtp/qmdatabase.h>
#include <votca/xtp/topology.h>
#include <boost/interprocess/sync/file_lock.hpp>

namespace votca { namespace xtp {



class StateSaverSQLite
{
public:
    StateSaverSQLite() { };
   ~StateSaverSQLite() { _db.Close(); }

    void Open(xtp::Topology &qmtop, const std::string &file, bool lock = true);
    void Close() { _db.Close(); }
    bool NextFrame();

    void WriteFrame();
    void WriteMeta(bool update);
    void WriteMolecules(bool update);
    void WriteSegTypes(bool update);
    void WriteSegments(bool update);
    void WriteFragments(bool update);
    void WriteAtoms(bool update);
    void WritePairs(bool update);
    void WriteSuperExchange(bool update);

    void ReadFrame();
    void ReadMeta(int topId);
    void ReadMolecules(int topId);
    void ReadSegTypes(int topId);
    void ReadSegments(int topId);
    void ReadFragments(int topId);
    void ReadAtoms(int topId);
    void ReadPairs(int topId);
    void ReadSuperExchange(int topId);

    int  FramesInDatabase();
    xtp::Topology *getTopology() { return _qmtop; }
    bool HasTopology(xtp::Topology *top);
    
    void LockStateFile();
    void UnlockStateFile();
    
private:
    xtp::Topology       *_qmtop;
    QMDatabase      _db;

    int             _frame;
    int             _current_frame;
    std::vector<int>     _frames;
    std::vector<int>     _topIds;

    std::string          _sqlfile;
    bool            _was_read;
    
    boost::interprocess::file_lock *_flock;
};

}}

#endif	// VOTCA_MD2QM_STATE_SAVER_SQLITE_H

