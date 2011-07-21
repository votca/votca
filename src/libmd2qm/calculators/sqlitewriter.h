#ifndef __VOTCA_MD2QM_SQLITEWRITE_H
#define __VOTCA_MD2QM_SQLITEWRITE_H

#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/statesaversqlite.h>

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

#endif

