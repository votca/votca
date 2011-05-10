#ifndef __VOTCA_MD2QM_SQLITEWRITE_H
#define __VOTCA_MD2QM_SQLITEWRITE_H

#include "qmcalculator.h"
#include "statesaversqlite.h"

class SQLiteWriter : public QMCalculator {
public:
    SQLiteWriter() {};
    virtual ~SQLiteWriter() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

private:
    StateSaverSQLite _saver;
};

inline void SQLiteWriter::Initialize(QMTopology *top, Property *options)
{
    _saver.Open(*top, "state.db");
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

