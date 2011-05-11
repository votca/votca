#include <statesaversqlite.h>

StateSaverSQLite::StateSaverSQLite()
    : _db(NULL)
{}

StateSaverSQLite::~StateSaverSQLite()
{
    Close();
}

void StateSaverSQLite::Open(QMTopology& qmtop, const string& file)
{
    int ret = sqlite3_open_v2(file.c_str(),&_db,SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("cannot open database");
    _frame = 0;
    _qmtop=&qmtop;
}

void StateSaverSQLite::Close()
{
    if(_db)
        sqlite3_close(_db);
}

void StateSaverSQLite::WriteFrame()
{
    sqlite3_stmt *stmt;
    char *error;

    int ret = sqlite3_exec(_db,
        "CREATE TABLE frames (_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "time REAL NOT NULL,"
        "step INT NOT NULL,"
        "box11 REAL NOT NULL,"
        "box12 REAL NOT NULL,"
        "box13 REAL NOT NULL,"
        "box21 REAL NOT NULL,"
        "box22 REAL NOT NULL,"
        "box23 REAL NOT NULL,"
        "box31 REAL NOT NULL,"
        "box32 REAL NOT NULL,"
        "box33 REAL NOT NULL)",
        NULL, NULL,  &error);
    if(ret != SQLITE_OK)
        throw std::runtime_error(string("cannot create frame table:\n") + error);


    ret = sqlite3_prepare_v2(_db,
            "INSERT INTO frames (time, step, box11, box12, box13, box21, box22, box23, box31, box32, box33) VALUES (?,?,?,?,?,?,?,?,?,?,?)"
            , -1, &stmt, NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("prepare insert frame statement failed");

    ret = sqlite3_bind_double(stmt, 1, _qmtop->getTime());
    if(ret != SQLITE_OK)
       throw std::runtime_error("bind failed 0 \n" + lexical_cast<string>(ret));

    sqlite3_bind_int(stmt, 2, _qmtop->getStep());
    for(int i=0; i<3; ++i)
        for(int j=0; j<9; ++j)
        sqlite3_bind_double(stmt, 3+i*3+j, _qmtop->getBox().get(i,j));

    sqlite3_step(stmt);
    ret = sqlite3_finalize(stmt);
    if(ret != SQLITE_OK)
       throw std::runtime_error("finalize failed\n");

    WriteMolecules();
    WriteCrgUnits();
    WriteBeads();
    WritePairs();

    _frame++;
}

void StateSaverSQLite::WriteMolecules()
{
    sqlite3_stmt *stmt;
    int ret = sqlite3_exec(_db,
        "CREATE TABLE molecules (_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL)",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("cannot create frame table:\n");


        ret = sqlite3_exec(_db,
        "BEGIN;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("begin transaction failed");

        
    ret = sqlite3_prepare_v2(_db,
            "INSERT INTO molecules (id, name) VALUES (?, ?)"
            , -1, &stmt, NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("prepare insert frame statement failed");


     int imol=0;
    for (MoleculeContainer::iterator iter = _qmtop->Molecules().begin();
            iter != _qmtop->Molecules().end() ; ++iter) {
        Molecule *mol=*iter;
        sqlite3_bind_int(stmt, 1, mol->getId());
        sqlite3_bind_text(stmt, 2, mol->getName().c_str(), -1, NULL);;
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }
    ret = sqlite3_exec(_db,
        "END;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("end transaction failed");


    ret = sqlite3_finalize(stmt);
    if(ret != SQLITE_OK)
       throw std::runtime_error("finalize failed\n");
}

void StateSaverSQLite::WriteCrgUnits() {
    sqlite3_stmt *stmt;
    int ret = sqlite3_exec(_db,
        "CREATE TABLE crgunits (_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "energy REAL NOT NULL,"
        "occ REAL NOT NULL)",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("cannot create frame table:\n");


        ret = sqlite3_exec(_db,
        "BEGIN;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("begin transaction failed");


    ret = sqlite3_prepare_v2(_db,
            "INSERT INTO crgunits (id, name, energy, occ) VALUES (?,?,?,?)"
            , -1, &stmt, NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("prepare insert frame statement failed");


     int imol=0;
    for (vector < QMCrgUnit *>::iterator iter = _qmtop->CrgUnits().begin(); iter!=_qmtop->CrgUnits().end(); ++iter) {
        QMCrgUnit *crg = *iter;
        sqlite3_bind_int(stmt, 1, crg->getId());
        sqlite3_bind_text(stmt, 2, crg->getName().c_str(), -1, NULL);;
        sqlite3_bind_double(stmt, 3, crg->getEnergy());
        sqlite3_bind_double(stmt, 4, crg->getOccupationProbability());
        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }
    ret = sqlite3_exec(_db,
        "END;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("end transaction failed");


    ret = sqlite3_finalize(stmt);
    if(ret != SQLITE_OK)
       throw std::runtime_error("finalize failed\n");
}

void StateSaverSQLite::WriteBeads() {
    sqlite3_stmt *stmt;
    int ret = sqlite3_exec(_db,
        "CREATE TABLE beads (_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "id INT NOT NULL,"
        "name TEXT NOT NULL,"
        "symmetry INT NOT NULL,"
        "type TEXT NOT NULL,"
        "resnr INT NOT NULL,"
        "mass REAL NOT NULL,"
        "charge REAL NOT NULL,"
        "crgunit INT NOT NULL,"
        "crgunit_index INT NOT NULL,"
        "pos_x REAL NOT NULL,"
        "pos_y REAL NOT NULL,"
        "pos_z REAL NOT NULL,"
        "u_x REAL NOT NULL,"
        "u_y REAL NOT NULL,"
        "u_z REAL NOT NULL,"
        "v_x REAL NOT NULL,"
        "v_y REAL NOT NULL,"
        "v_z REAL NOT NULL,"
        "molid INT NOT NULL)",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("cannot create frame table:\n");


        ret = sqlite3_exec(_db,
        "BEGIN;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("begin transaction failed");


    ret = sqlite3_prepare_v2(_db,
            "INSERT INTO beads (id,name,symmetry,type,resnr,mass,charge,crgunit,"
            "crgunit_index,pos_x,pos_y,pos_z,u_x,u_y,u_z,v_x,v_y,v_z,molid) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
            , -1, &stmt, NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("prepare insert frame statement failed");

    for (BeadContainer::iterator iter = _qmtop->Beads().begin();
            iter != _qmtop->Beads().end(); ++iter) {
        QMBead *bi = dynamic_cast<QMBead*> (*iter);
        sqlite3_bind_int(stmt,     1, bi->getId());
        sqlite3_bind_text(stmt,    2, bi->getName().c_str(), -1, NULL);;
        sqlite3_bind_int(stmt,     3,(int)bi->getSymmetry());
        sqlite3_bind_text(stmt,    4, bi->getType()->getName().c_str(), -1, NULL);;
        sqlite3_bind_int(stmt,     5, bi->getResnr());
        sqlite3_bind_double(stmt,  6, bi->getM());
        sqlite3_bind_double(stmt,  7, bi->getQ());
        sqlite3_bind_int(stmt,     8, bi->GetCrgUnit()->getId());;
        sqlite3_bind_int(stmt,     9, bi->getiPos());
        
        sqlite3_bind_double(stmt, 10, bi->getPos().getX());
        sqlite3_bind_double(stmt, 11, bi->getPos().getY());
        sqlite3_bind_double(stmt, 12, bi->getPos().getZ());
        sqlite3_bind_double(stmt, 13, bi->getU().getX());
        sqlite3_bind_double(stmt, 14, bi->getU().getY());
        sqlite3_bind_double(stmt, 15, bi->getU().getZ());
        sqlite3_bind_double(stmt, 16, bi->getV().getX());
        sqlite3_bind_double(stmt, 17, bi->getV().getY());
        sqlite3_bind_double(stmt, 18, bi->getV().getZ());
        
        sqlite3_bind_int(stmt,    19, bi->getMolecule()->getId());

        sqlite3_step(stmt);
        sqlite3_reset(stmt);


    }
        ret = sqlite3_exec(_db,
        "END;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("end transaction failed");


    ret = sqlite3_finalize(stmt);
    if(ret != SQLITE_OK)
       throw std::runtime_error("finalize failed\n");

}

void StateSaverSQLite::WritePairs() {
    sqlite3_stmt *stmt;
    int ret = sqlite3_exec(_db,
        "CREATE TABLE pairs (_id INTEGER PRIMARY KEY AUTOINCREMENT,"
        "crgunit1 INT NOT NULL,"
        "crgunit2 INT NOT NULL,"
        "rate12 REAL NOT NULL,"
        "rate21 REAL NOT NULL,"
        "r_x REAL NOT NULL,"
        "r_y REAL NOT NULL,"
        "r_z REAL NOT NULL)",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("cannot create pairs table:\n");


        ret = sqlite3_exec(_db,
        "BEGIN;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("begin transaction failed");


    ret = sqlite3_prepare_v2(_db,
            "INSERT INTO pairs (crgunit1, crgunit2, rate12, rate21, r_x, r_y,r_z)"
            " VALUES (?,?,?,?,?,?,?)"
            , -1, &stmt, NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("prepare insert frame statement failed");

    QMNBList &nblist = _qmtop->nblist();
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        QMPair *pair = *iter;
        QMCrgUnit *crg1 = (*iter)->first;
        QMCrgUnit *crg2 = (*iter)->second;

        sqlite3_bind_int(stmt,     1, (int)crg1->getId());
        sqlite3_bind_int(stmt,     2, (int)crg2->getId());
        sqlite3_bind_double(stmt,  3, pair->rate12());
        sqlite3_bind_double(stmt,  4, pair->rate21());
        sqlite3_bind_double(stmt,  5, pair->r().getX());
        sqlite3_bind_double(stmt,  6, pair->r().getY());
        sqlite3_bind_double(stmt,  7, pair->r().getZ());

        sqlite3_step(stmt);
        sqlite3_reset(stmt);
    }

    ret = sqlite3_exec(_db,
        "END;",
        NULL, NULL,  NULL);
    if(ret != SQLITE_OK)
        throw std::runtime_error("end transaction failed");


    ret = sqlite3_finalize(stmt);
    if(ret != SQLITE_OK)
       throw std::runtime_error("finalize failed\n");

}
