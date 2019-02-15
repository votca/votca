/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#include <votca/xtp/qmdatabase.h>

namespace votca {
namespace xtp {

void QMDatabase::onCreate() {
  // Table format frames
  Exec(
      "CREATE TABLE frames ("
      "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
      "id       INT NOT NULL,"
      "time     REAL NOT NULL,"
      "step     INT NOT NULL,"
      "box11    REAL NOT NULL,"
      "box12    REAL NOT NULL,"
      "box13    REAL NOT NULL,"
      "box21    REAL NOT NULL,"
      "box22    REAL NOT NULL,"
      "box23    REAL NOT NULL,"
      "box31    REAL NOT NULL,"
      "box32    REAL NOT NULL,"
      "box33    REAL NOT NULL,"
      "canRigid INT NOT NULL)");

  // Table format molecules
  Exec(
      "CREATE TABLE molecules ("
      "_id    INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame  INT NOT NULL,"
      "top    INT NOT NULL,"
      "id     INT NOT NULL,"
      "name   TEXT NOT NULL,"
      "type   TEXT NOT NULL)");

  // Table format segments
  Exec(
      "CREATE TABLE segments ("
      "_id        INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame  INT NOT NULL,"
      "top    INT NOT NULL,"
      "id     INT NOT NULL,"
      "name   TEXT NOT NULL,"
      "type   TEXT NOT NULL,"
      "mol    INT NOT NULL,"

      "posX   REAL NOT NULL,"
      "posY   REAL NOT NULL,"
      "posZ   REAL NOT NULL,"

      "UnCnNe     REAL DEFAULT 0,"
      "UnCnNh     REAL DEFAULT 0,"
      "UcNcCe     REAL DEFAULT 0,"
      "UcNcCh     REAL DEFAULT 0,"
      "UcCnNe     REAL DEFAULT 0,"
      "UcCnNh     REAL DEFAULT 0,"
      "UnXnNs     REAL DEFAULT 0,"
      "UnXnNt     REAL DEFAULT 0,"
      "UxNxXs     REAL DEFAULT 0,"
      "UxNxXt     REAL DEFAULT 0,"
      "UxXnNs     REAL DEFAULT 0,"
      "UxXnNt     REAL DEFAULT 0,"
      "eAnion     REAL DEFAULT 0,"
      "eNeutral   REAL DEFAULT 0,"
      "eCation    REAL DEFAULT 0,"
      "eSinglet     REAL DEFAULT 0,"
      "eTriplet     REAL DEFAULT 0,"

      "has_e      INT  DEFAULT 0,"
      "has_h      INT  DEFAULT 0,"
      "has_s      INT  DEFAULT 0,"
      "has_t      INT  DEFAULT 0,"

      "occPe      REAL DEFAULT -1,"
      "occPh      REAL DEFAULT -1,"
      "occPs      REAL DEFAULT -1,"
      "occPt      REAL DEFAULT -1)");

  // Table format segment types
  Exec(
      "CREATE TABLE segmentTypes ("
      "_id       INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame     INT NOT NULL,"
      "top       INT NOT NULL,"
      "id        INT NOT NULL,"
      "name      TEXT NOT NULL,"
      "basis     TEXT NOT NULL,"
      "orbfile   TEXT NOT NULL,"
      "torbnrs   TEXT NOT NULL,"
      "coordfile TEXT NOT NULL,"
      "canRigid  INT NOT NULL)");

  // Table format fragments
  Exec(
      "CREATE TABLE fragments ("
      "_id    INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame  INT NOT NULL,"
      "top    INT NOT NULL,"
      "id     INT NOT NULL,"
      "name   TEXT NOT NULL,"
      "type   TEXT NOT NULL,"
      "mol    INT NOT NULL,"
      "seg    INT NOT NULL,"

      "posX   REAL NOT NULL,"
      "posY   REAL NOT NULL,"
      "posZ   REAL NOT NULL,"
      "symmetry INT NOT NULL,"
      "leg1   INT NOT NULL,"
      "leg2   INT NOT NULL,"
      "leg3   INT NOT NULL)");

  // Table format atoms
  Exec(
      "CREATE TABLE atoms ("
      "_id   INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame INT NOT NULL,"
      "top   INT NOT NULL,"
      "id    INT NOT NULL,"
      "name  TEXT NOT NULL,"
      "type  INT NOT NULL,"

      "mol   INT NOT NULL,"
      "seg   INT NOT NULL,"
      "frag  INT NOT NULL,"

      "resnr   INT NOT NULL,"
      "resname TEXT NOT NULL,"

      "posX    REAL NOT NULL,"
      "posY    REAL NOT NULL,"
      "posZ    REAL NOT NULL,"
      "weight  REAL NOT NULL,"
      "element TEXT NOT NULL,"
      "qmid    INT NOT NULL,"
      "qmPosX  REAL NOT NULL,"
      "qmPosY  REAL NOT NULL,"
      "qmPosZ  REAL NOT NULL)");

  // Table format pairs
  Exec(
      "CREATE TABLE pairs ("
      "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame      INT NOT NULL,"
      "top        INT NOT NULL,"
      "id         INT NOT NULL,"

      "seg1       INT NOT NULL,"
      "seg2       INT NOT NULL,"

      "drX        REAL NOT NULL,"
      "drY        REAL NOT NULL,"
      "drZ        REAL NOT NULL,"

      "lOe        REAL DEFAULT 0,"
      "lOh        REAL DEFAULT 0,"
      "lOs        REAL DEFAULT 0,"
      "lOt        REAL DEFAULT 0,"

      "has_e      INT  DEFAULT 0,"
      "has_h      INT  DEFAULT 0,"
      "has_s      INT  DEFAULT 0,"
      "has_t      INT  DEFAULT 0,"
      "rate12e    REAL DEFAULT 0,"
      "rate21e    REAL DEFAULT 0,"
      "rate12h    REAL DEFAULT 0,"
      "rate21h    REAL DEFAULT 0,"
      "rate12s    REAL DEFAULT 0,"
      "rate21s    REAL DEFAULT 0,"
      "rate12t    REAL DEFAULT 0,"
      "rate21t    REAL DEFAULT 0,"
      "Jeff2e     REAL DEFAULT 0,"
      "Jeff2h     REAL DEFAULT 0,"
      "Jeff2s     REAL DEFAULT 0,"
      "Jeff2t     REAL DEFAULT 0,"
      "type       INT  DEFAULT 0)");

  // Super-exchange types
  Exec(
      "CREATE TABLE superExchange ("
      "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
      "frame      INT NOT NULL,"
      "top        INT NOT NULL,"
      "type      TEXT NOT NULL)");
}

}  // namespace xtp
}  // namespace votca
