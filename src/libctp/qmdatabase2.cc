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
        "id     INT NOT NULL,"
        "time   REAL NOT NULL,"
        "step   INT NOT NULL,"
        "box11  REAL NOT NULL,"
        "box12  REAL NOT NULL,"
        "box13  REAL NOT NULL,"
        "box21  REAL NOT NULL,"
        "box22  REAL NOT NULL,"
        "box23  REAL NOT NULL,"
        "box31  REAL NOT NULL,"
        "box32  REAL NOT NULL,"
        "box33  REAL NOT NULL)");
   
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

        "enerA  REAL,"
        "enerN  REAL,"
        "enerC  REAL,"
        "occ    REAL)");

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
        "symmetry INT NOT NULL)");

    // Table format atoms
    Exec("CREATE TABLE atoms ("
         "_id   INTEGER PRIMARY KEY AUTOINCREMENT,"
         "frame INT NOT NULL,"
         "top   INT NOT NULL,"
         "id    INT NOT NULL,"
         "name  TEXT NOT NULL,"
         "type  TEXT NOT NULL,"

         "mol   INT NOT NULL,"
         "seg   INT NOT NULL,"
         "frag  INT NOT NULL,"

         "resnr   INT NOT NULL,"
         "resname TEXT NOT NULL,"

         "posX    REAL NOT NULL,"
         "posY    REAL NOT NULL,"
         "posZ    REAL NOT NULL,"
         "weight  REAL NOT NULL)");
               
    // Table format pairs
    Exec("CREATE TABLE pairs ("
        "_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "conjseg1 INT NOT NULL,"
        "conjseg2 INT NOT NULL,"
        "rate12 REAL NOT NULL,"
        "rate21 REAL NOT NULL,"
        "r_x REAL NOT NULL,"
        "r_y REAL NOT NULL,"
        "r_z REAL NOT NULL,"
        "deleted INT DEFAULT 0)");

    // Additional pair properties
    Exec("CREATE TABLE pair_properties ("
        "_id    INTEGER PRIMARY KEY AUTOINCREMENT,"
        "pair   INTEGER NOT NULL,"
        "key    TEXT NOT NULL,"
        "value  REAL NOT NULL)");

    // Table format coupling elements
    Exec("CREATE TABLE pair_integrals ("
        "pair   INTEGER NOT NULL,"
        "num    INTEGER NOT NULL,"
        "J      REAL NOT NULL)");
}

}}
