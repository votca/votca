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

#include <votca/ctp/qmdatabase2.h>

namespace votca { namespace ctp {

void QMDatabase2::onCreate()
{
    // Table format frames
    Exec("CREATE TABLE frames ("
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
    Exec("CREATE TABLE molecules ("
        "_id    INTEGER PRIMARY KEY AUTOINCREMENT,"
        "frame  INT NOT NULL,"
        "top    INT NOT NULL,"
        "id     INT NOT NULL,"
        "name   TEXT NOT NULL,"
        "type   TEXT NOT NULL)");
      
    // Table format segments
    Exec("CREATE TABLE segments ("
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

        "lI_AN      REAL DEFAULT 0,"
        "lI_NA      REAL DEFAULT 0,"
        "lI_CN      REAL DEFAULT 0,"
        "lI_NC      REAL DEFAULT 0,"
        "eI_A       REAL DEFAULT 0,"
        "eI_C       REAL DEFAULT 0,"
        "eAnion     REAL DEFAULT 0,"
        "eNeutral   REAL DEFAULT 0,"
        "eCation    REAL DEFAULT 0,"

        "occPe      REAL DEFAULT -1,"
        "occPh      REAL DEFAULT -1)");

    // Table format segment types
    Exec("CREATE TABLE segmentTypes ("
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
    Exec("CREATE TABLE fragments ("
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
    Exec("CREATE TABLE atoms ("
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
    Exec("CREATE TABLE pairs ("
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

        "has_e      INT  DEFAULT 0,"
        "has_h      INT  DEFAULT 0,"
        "rate12e    REAL DEFAULT 0,"
        "rate21e    REAL DEFAULT 0,"
        "rate12h    REAL DEFAULT 0,"
        "rate21h    REAL DEFAULT 0,"
        "Jeff2e     REAL DEFAULT 0,"
        "Jeff2h     REAL DEFAULT 0)");


}

}}
