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

#ifndef __VOTCA_MD2QM_SQLITEWRITE_H
#define __VOTCA_MD2QM_SQLITEWRITE_H

#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/statesaversqlite.h>

namespace votca { namespace ctp {

/** \brief Writes out data to sqlite3 format

Callname: sqlitewriter

Writes out the following pair properties to sqlite3 format:
 - Distance between sites i and j in the pair in nm.
 - Distance vector coordinates x,y,z in nm.
 - Transfer integral in eV (or the effective value if more frontier orbitals are involved)
 - Intramolecular reorganization energy in eV.
 - Outer sphere reorganization energy in eV.

Writes out the following site properties to sqlite3 format:
 - Site energy in eV
 - Site occupation probability
 - Coordinates x,y,z of the center of mass in nm
*/
class SQLiteWriter : public QMCalculator {
public:
    SQLiteWriter() {};
    virtual ~SQLiteWriter() {};

    const char *Description() { return "Writes out data to sqlite3 format"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

private:
    StateSaverSQLite _saver;
};

inline void SQLiteWriter::Initialize(QMTopology *top, Property *options)
{
    _saver.Open(*top, "state_calculator.db");
}

inline bool SQLiteWriter::EvaluateFrame(QMTopology *top)
{
    _saver.WriteFrame();
    return true;
}

inline void SQLiteWriter::EndEvaluate(QMTopology *top)
{
    _saver.Close();
}

}}

#endif

