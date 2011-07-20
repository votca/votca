/* 
 * File:   StateSaverSQLite.h
 * Author: ruehle
 *
 * Created on April 11, 2011, 2:00 PM
 */

#ifndef __VOTCA_MD2QM_STATESAVERSQLITE_H
#define	__VOTCA_MD2QM_STATESAVERSQLITE_H

#include <stdio.h>
#include <map>
#include "qmbead.h"
#include "qmtopology.h"
#include "qmpair.h"
#include "qmnblist.h"
#include "qmdatabase.h"

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

    int FramesInDatabase() { return _frames.size(); }
    
private:
    int _frame;
    QMDatabase _db;
    
    void WriteMolecules(int frameid);
    void WriteConjugatedSegments(int frameid);
    void WriteBeads(int frameid);
    void WritePairs(int frameid);
    void WriteIntegrals(QMPair *pair);

    void ReadFrame(void);
    void ReadMolecules(void);
    void ReadConjugatedSegments(void);
    void ReadBeads(void);
    void ReadPairs(void);
    void ReadIntegrals();

    template<typename T>
    void WriteCustomProperties(string table, std::map<string, T> &properties,
        const string field_key, const string field_value);
    
    QMTopology *_qmtop;

    vector<int> _frames;
    int _current_frame;

    map<int,int> _conjseg_id_map;
    bool _was_read;
};

#endif	/* __VOTCA_MD2QM_StateSaverSQLite_H */

