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
#include <votca/tools/database.h>

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
private:
    int _frame;
    Database _db;
    
    void WriteMolecules(int frameid);
    void WriteCrgUnits(int frameid);
    void WriteBeads(int frameid);
    void WritePairs(int frameid);

    void ReadFrame(void);
    void ReadMolecules(void);
    void ReadCrgUnits(void);
    void ReadBeads(void);
    void ReadPairs(void);

    QMTopology *_qmtop;

    vector<int> _frames;
    int _current_frame;
};

#endif	/* __VOTCA_MD2QM_StateSaverSQLite_H */

