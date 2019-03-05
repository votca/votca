/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#ifndef VOTCA_XTP_STATE_SAVER_H
#define VOTCA_XTP_STATE_SAVER_H

#include <boost/interprocess/sync/file_lock.hpp>
#include <map>
#include <stdio.h>
#include <votca/xtp/qmdatabase.h>
#include <votca/xtp/topology.h>

namespace votca {
namespace xtp {

class StateSaverSQLite {
 public:
  StateSaverSQLite(){};
  ~StateSaverSQLite() { _db.Close(); }

  void Open(const std::string &file, bool lock = true);
  void Close() { _db.Close(); }
  bool NextFrame();

  void WriteFrame(const Topology &top);
  void WriteMeta(const Topology &top, bool update);
  void WriteSegments(const Topology &top, bool update);
  void WriteAtoms(const Topology &top, bool update);
  void WritePairs(const Topology &top, bool update);

  void ReadFrame(Topology &top);
  void ReadMeta(Topology &top, int topId);
  void ReadSegments(Topology &top, int topId);
  void ReadAtoms(Topology &top, int topId);
  void ReadPairs(Topology &top, int topId);

  int FramesInDatabase();
  void LockStateFile();
  void UnlockStateFile();

 private:
  bool HasTopology(const Topology &top);
  QMDatabase _db;

  int _frame;
  int _current_frame;
  std::vector<int> _frames;
  std::vector<int> _topIds;

  std::string _sqlfile;
  bool _was_read;

  boost::interprocess::file_lock _flock;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_MD2QM_STATE_SAVER_SQLITE_H
