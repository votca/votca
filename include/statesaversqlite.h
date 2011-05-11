/* 
 * File:   StateSaverSQLite.h
 * Author: ruehle
 *
 * Created on April 11, 2011, 2:00 PM
 */

#ifndef __VOTCA_MD2QM_STATESAVERSQLITE_H
#define	__VOTCA_MD2QM_STATESAVERSQLITE_H

#include <stdio.h>
#include "qmbead.h"
#include "qmtopology.h"
#include "qmpair.h"
#include "qmnblist.h"
#include <sqlite3.h>

class StateSaverSQLite
{
public:
    StateSaverSQLite();
    ~StateSaverSQLite();

    void Open(QMTopology & qmtop,const string &file);
    void Close();
    void WriteFrame();
private:
    int _frame;
    sqlite3 *_db;
    
    void WriteMolecules();
    void WriteCrgUnits();
    void WriteBeads();
    void WritePairs();

    QMTopology *_qmtop;
};

#endif	/* __VOTCA_MD2QM_StateSaverSQLite_H */

