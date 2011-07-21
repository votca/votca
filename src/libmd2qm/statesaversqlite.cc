#include <statesaversqlite.h>
#include <votca/tools/statement.h>

StateSaverSQLite::StateSaverSQLite()
{}

StateSaverSQLite::~StateSaverSQLite()
{
    Close();
}

void StateSaverSQLite::Open(QMTopology& qmtop, const string& file)
{
    _db.OpenHelper(file.c_str());
    _frame = 0;
    _qmtop=&qmtop;

    // now query available frames in database
    Statement *stmt = _db.Prepare("SELECT _id FROM frames;");
    while(stmt->Step() != SQLITE_DONE) {
        _frames.push_back(stmt->Column<int>(0));
    }
    delete stmt;

    cout << "Found " << _frames.size() << " in database\n";
    _current_frame = -1;
    _was_read=false;
}

void StateSaverSQLite::Close()
{
    _db.Close();
}

bool StateSaverSQLite::NextFrame()
{
    _qmtop->Cleanup();
    _qmtop->nblist().Cleanup();
    _qmtop->CreateResidue("dummy");

    _current_frame++;
    if(_current_frame >= _frames.size())
        return false;

    ReadFrame();
    ReadMolecules();
    ReadConjugatedSegments();
    ReadBeads();
    ReadPairs();
    _was_read=true;
    return true;
}

void StateSaverSQLite::ReadFrame(void)
{
    cout << "reading frame " << _current_frame+1 << ", id " << _frames[_current_frame] << endl;
    Statement *stmt =
        _db.Prepare("SELECT time, step, "
            "box11, box12, box13, "
            "box21, box22, box23, "
            "box31, box32, box33 "
            "FROM frames WHERE _id = ?;");
    stmt->Bind(1, _frames[_current_frame]);

    if(stmt->Step() == SQLITE_DONE)
        throw std::runtime_error("cannot read frame, database might have changed while program is running");

    _qmtop->setTime(stmt->Column<double>(0));
    _qmtop->setStep(stmt->Column<int>(1));
    cout << "time " << stmt->Column<double>(0) << endl;
    matrix m;
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            m.set(i, j, stmt->Column<double>(2+i*3+j));
    _qmtop->setBox(m);
    _qmtop->setDatabaseId(_frames[_current_frame]);
    delete stmt;
}

void StateSaverSQLite::WriteFrame()
{
    _db.BeginTransaction();
    Statement *stmt;
    
    // create a new entry or update existing one
    if(_qmtop->getDatabaseId() == 0) {
        stmt = _db.Prepare(
            "INSERT INTO frames (time, step, "
            "box11, box12, box13, box21, box22, box23, box31, box32, box33) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?)");
    } else {
        stmt = _db.Prepare(
            "UPDATE frames SET time = ?, step = ?, "
            "box11 = ?, box12 = ?, box13 = ?, box21 = ?, box22 = ?, box23 = ?, box31 = ?, box32 = ?, box33 = ? "
            "WHERE _id = ?");
        stmt->Bind(12, _qmtop->getDatabaseId());
    }

    stmt->Bind(1, _qmtop->getTime());
    stmt->Bind(2, _qmtop->getStep());
    matrix m = _qmtop->getBox();
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            stmt->Bind(3+i*3+j, _qmtop->getBox().get(i,j));

    stmt->Step();
    if(_qmtop->getDatabaseId() == 0)
        _qmtop->setDatabaseId(_db.LastInsertRowId());

    delete stmt;


    int frameid = _db.LastInsertRowId();
    WriteMolecules(_qmtop->getDatabaseId());
    WriteConjugatedSegments(_qmtop->getDatabaseId());
    WriteBeads(_qmtop->getDatabaseId());
    WritePairs(_qmtop->getDatabaseId());

    _db.EndTransaction();
    _frame++;
    _conjseg_id_map.clear();
}

void StateSaverSQLite::WriteMolecules(int frameid)
{
    Statement *stmt;
    if(_was_read) return;
    stmt = _db.Prepare("INSERT INTO molecules (frame, id, name) VALUES (?, ?, ?)");

    int imol=0;
    stmt->Bind(1, frameid);

    for (MoleculeContainer::iterator iter = _qmtop->Molecules().begin();
            iter != _qmtop->Molecules().end() ; ++iter) {
        Molecule *mol=*iter;
        stmt->Bind(2, mol->getId());
        stmt->Bind(3, mol->getName());
        stmt->Step();
        stmt->Reset();
        //mol->setDBId(_db.LastInsertRowId());
    }

    delete stmt;
}

void StateSaverSQLite::ReadMolecules(void)
{
   Statement *stmt =
        _db.Prepare("SELECT name FROM molecules WHERE frame = ?;");
    stmt->Bind(1, _frames[_current_frame]);

    while(stmt->Step() != SQLITE_DONE) {
        Molecule *mol = _qmtop->CreateMolecule(stmt->Column<string>(0));
    }
    delete stmt;
}

void StateSaverSQLite::WriteConjugatedSegments(int frameid) {
    Statement *update_stmt = _db.Prepare(
            "UPDATE conjsegs SET name = ?, type = ?, molecule = ?, frame = ? WHERE _id = ?"
    );
    Statement *insert_stmt = _db.Prepare(
            "INSERT INTO conjsegs (name, type, molecule, frame) VALUES (?,?,?,?)"
    );


    int imol=0;
    for (vector < QMCrgUnit *>::iterator iter = _qmtop->CrgUnits().begin(); iter!=_qmtop->CrgUnits().end(); ++iter) {
        QMCrgUnit *crg = *iter;
        Statement *stmt;
        if(crg->getInDatabase()) {
            stmt = update_stmt;
            stmt->Bind(5, (int)crg->getId());
        }
        else
            stmt = insert_stmt;

        // cout << "conjseg " << crg->getId() << " " << crg->GetCom() << endl;

        stmt->Bind(1, crg->getName());;
        stmt->Bind(2, crg->getType()->GetName());;
        stmt->Bind<int>(3, crg->getMolId());;
        stmt->Bind(4, frameid);;
        stmt->Step();
        stmt->Reset();

        if(crg->getInDatabase())
            _conjseg_id_map[crg->getId()] = crg->getId();
        else
            _conjseg_id_map[crg->getId()] = _db.LastInsertRowId();

        WriteCustomProperties(_conjseg_id_map[crg->getId()], crg->DoubleValues(), "conjseg_properties", "conjseg");
    }
    delete insert_stmt;
    delete update_stmt;
}

void StateSaverSQLite::ReadConjugatedSegments(void)
{
   Statement *stmt =
        _db.Prepare("SELECT _id, name, type, molecule FROM conjsegs WHERE frame = ?;");
    stmt->Bind(1, _frames[_current_frame]);

    while(stmt->Step() != SQLITE_DONE) {
        QMCrgUnit *acrg = _qmtop->CreateCrgUnit(stmt->Column<int>(0),
                stmt->Column<string>(1), stmt->Column<string>(2), stmt->Column<int>(3));
        acrg->setInDatabase(true);
        ReadCustomProperties(acrg->getId(), acrg->DoubleValues(), "conjseg_properties", "conjseg");
    }
    delete stmt;
}

void StateSaverSQLite::WriteBeads(int frameid) {
    if(_was_read) return;
    Statement *stmt= _db.Prepare(
            "INSERT INTO rigidfrags (id,name,symmetry,type,resnr,mass,charge,conjseg_id,"
            "conjseg_index,pos_x,pos_y,pos_z,u_x,u_y,u_z,v_x,v_y,v_z,molecule,frame) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");

    int i = 0;
    for (BeadContainer::iterator iter = _qmtop->Beads().begin();
            iter != _qmtop->Beads().end(); ++iter) {
        QMBead *bi = dynamic_cast<QMBead*> (*iter);
        stmt->Bind(1, bi->getId());
        stmt->Bind(2, bi->getName());
        stmt->Bind(3,(int)bi->getSymmetry());
        stmt->Bind(4, bi->getType()->getName());
        stmt->Bind(5, bi->getResnr());
        stmt->Bind(6, bi->getM());
        stmt->Bind(7, bi->getQ());
        stmt->Bind<int>(8, _conjseg_id_map[bi->GetCrgUnit()->getId()]);
        stmt->Bind(9, bi->getiPos());
        
        stmt->Bind(10, bi->getPos().getX());
        stmt->Bind(11, bi->getPos().getY());
        stmt->Bind(12, bi->getPos().getZ());
        stmt->Bind(13, bi->getU().getX());
        stmt->Bind(14, bi->getU().getY());
        stmt->Bind(15, bi->getU().getZ());
        stmt->Bind(16, bi->getV().getX());
        stmt->Bind(17, bi->getV().getY());
        stmt->Bind(18, bi->getV().getZ());
        stmt->Bind(19, bi->getMolecule()->getId());
        stmt->Bind(20, frameid);
        
        stmt->Step();
        stmt->Reset();
    }
}

void StateSaverSQLite::ReadBeads() {
    Statement *stmt =
        _db.Prepare("SELECT "
            "rigidfrags._id, rigidfrags.id, rigidfrags.name, symmetry, rigidfrags.type, resnr, mass, charge, molecule, "
            "conjseg_id, conjseg_index, pos_x, pos_y, pos_z, "
            "u_x, u_y, u_z, v_x, v_y, v_z "
            "FROM rigidfrags WHERE (frame = ?)");
    stmt->Bind(1, _frames[_current_frame]);

    while (stmt->Step() != SQLITE_DONE) {
        int id = stmt->Column<int>(1);
        string bead_name = stmt->Column<string>(2);
        byte_t symmetry = stmt->Column<int>(3);
        string type_name = stmt->Column<string>(4);

        int resnr = stmt->Column<int>(5);
        double M =  stmt->Column<double>(6);
        double Q =  stmt->Column<double>(7);
        int molecule = stmt->Column<double>(8);

        int conjseg_id =  stmt->Column<int>(9);
        int conjseg_index =  stmt->Column<int>(10);

 	vec pos(stmt->Column<double>(11), stmt->Column<double>(12), stmt->Column<double>(13));
        vec u(stmt->Column<double>(14), stmt->Column<double>(15), stmt->Column<double>(16));
        vec v(stmt->Column<double>(17), stmt->Column<double>(18), stmt->Column<double>(19));
        //cout << "read bead with conksegment " << conjseg_id << endl;

        BeadType *type = _qmtop->GetOrCreateBeadType(type_name);
	resnr = 0;
        QMBead *bead = dynamic_cast<QMBead*>(_qmtop->CreateBead(symmetry, bead_name, type, resnr, M, Q));
        _qmtop->getMolecule(molecule)->AddBead(bead, bead_name);

        QMCrgUnit * acrg = _qmtop->getCrgUnit(conjseg_id);
        if(acrg == NULL)
            throw std::runtime_error("error reading rigid fragments: charge unit not found");
        //acrg = _qmtop->CreateCrgUnit(crg_unit_name, type_name, molid);

        bead->setCrg(acrg);
        bead->setiPos(conjseg_index);
        bead->setPos(pos);
        bead->setU(u);
        bead->setV(v);
        //cout << "pos " << pos << endl;
        bead->UpdateCrg();
    }
    delete stmt;
}

void StateSaverSQLite::WritePairs(int frameid) {
    Statement *stmt;

    _db.Exec("UPDATE pairs SET deleted=1 WHERE conjseg1 IN (SELECT _id from conjsegs WHERE frame = " 
        + lexical_cast<string>(frameid) +  ")");
    //Statement *update_stmt = _db.Prepare(
    //        "UPDATE pairs SET rate12 = ?, rate21 = ?, r_x = ?, r_y = ?,r_z = ?, deleted = 0 WHERE conjseg1 = ? AND conjseg2 = ?");
    Statement *update_stmt = _db.Prepare(
            "UPDATE pairs SET rate12 = ?, rate21 = ?, r_x = ?, r_y = ?,r_z = ?, deleted = 0, conjseg1 = ?, conjseg2 = ? WHERE _id = ?");
    Statement *insert_stmt = _db.Prepare(
            "INSERT INTO pairs (rate12, rate21, r_x, r_y,r_z, conjseg1, conjseg2)"
            " VALUES (?,?,?,?,?,?,?)");

    int i=0;
QMNBList &nblist = _qmtop->nblist();
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        QMPair *pair = *iter;
        QMCrgUnit *crg1 = (*iter)->first;
        QMCrgUnit *crg2 = (*iter)->second;

        if(pair->getInDatabase())
            stmt = update_stmt;
        else
            stmt = insert_stmt;
        //cout << pair->getId() << " " << _conjseg_id_map[crg1->getId()] << " " << _conjseg_id_map[crg2->getId()] <<
        // " " << pair->rate12() << " " << pair->rate21() << " " <<  pair->r() << endl;
        stmt->Bind(1, pair->rate12());
        stmt->Bind(2, pair->rate21());
        stmt->Bind(3, pair->r().getX());
        stmt->Bind(4, pair->r().getY());
        stmt->Bind(5, pair->r().getZ());
        stmt->Bind(6, _conjseg_id_map[crg1->getId()]);
        stmt->Bind(7, _conjseg_id_map[crg2->getId()]);
        if(pair->getInDatabase())
            stmt->Bind(8, pair->getId());
        ++i;
        stmt->Step();
        if(!pair->getInDatabase())
            pair->setId(_db.LastInsertRowId());
        pair->setInDatabase(true);
        stmt->Reset();

        WriteIntegrals(pair);
        WriteCustomProperties(pair->getId(), pair->DoubleValues(), "pair_properties", "pair");
    }
    // cout << "written pairs: " << i << endl;
    _db.Exec("DELETE FROM pairs WHERE deleted=1");
    delete update_stmt;
    delete insert_stmt;
}

void StateSaverSQLite::ReadPairs(void)
{
    Statement *stmt =
    _db.Prepare("SELECT pairs._id, conjseg1, conjseg2, r_x, r_y, r_z FROM pairs,conjsegs "
            "WHERE (conjsegs.frame = ? and  conjseg1 = conjsegs._id)");
    stmt->Bind(1, _frames[_current_frame]);
    while (stmt->Step() != SQLITE_DONE) {
        int id1 = stmt->Column<int>(1);
        int id2 = stmt->Column<int>(2);
        QMCrgUnit *crg1 = _qmtop->getCrgUnit(id1);
        QMCrgUnit *crg2 = _qmtop->getCrgUnit(id2);
        if(crg1 == NULL)
            throw std::runtime_error("broken database, pair refers to non-existent conjugated segment");
        if(crg2 == NULL)
            throw std::runtime_error("broken database, pair refers to non-existent conjugated segment");

        QMPair *pair = new QMPair(crg1, crg2, _qmtop);
        _qmtop->nblist().AddPair(pair);
        vec r1(stmt->Column<double>(3), stmt->Column<double>(4), stmt->Column<double>(5));
        vec r2 = pair->r();
        if(abs(r2 - r1) > 1e-6)
            cerr << "WARNING: pair (" << id1 << ", " << id2 << ") distance differs by more than 1e-6 from the value in the database\nread: " << r1 << " calculated: " << r2 << endl ;
        pair->setInDatabase(true);
        pair->setId(stmt->Column<int>(0));
        ReadCustomProperties(pair->getId(), pair->DoubleValues(), "pair_properties", "pair");
        // cout << stmt->Column<int>(0) << "foo" << endl;
    }
    // cout << "read pairs: " + _qmtop->nblist().size() << endl;
    ReadIntegrals();
    delete stmt;
}

void StateSaverSQLite::WriteIntegrals(QMPair *pair)
{
    _db.Exec("DELETE FROM pair_integrals WHERE pair = " + lexical_cast<string>(pair->getId()));

    Statement *stmt =
    _db.Prepare("INSERT INTO pair_integrals (pair, num, J) VALUES (?, ?, ?)");

    for(int i=0; i<pair->Js().size();++i) {
        stmt->Bind(1, pair->getId());
        stmt->Bind(2, i);
        stmt->Bind(3, pair->Js()[i]);
        stmt->Step();
        stmt->Reset();
    }

    delete stmt;
}

void StateSaverSQLite::ReadIntegrals()
{
    Statement *stmt =
    _db.Prepare("SELECT J FROM pair_integrals WHERE pair = ? ORDER BY num");

    for(QMNBList::iterator iter=_qmtop->nblist().begin(); iter!=_qmtop->nblist().end(); ++iter) {
        stmt->Bind(1, (*iter)->getId());
        while (stmt->Step() != SQLITE_DONE)
            (*iter)->Js().push_back(stmt->Column<double>(0));
        stmt->Reset();
    }
    delete stmt;
}


template<typename T>
void StateSaverSQLite::WriteCustomProperties(int object_id, std::map<string, T> &properties,
        string table, const string field_objectid, const string field_key, const string field_value)
{
    _db.Exec("DELETE FROM " + table + " WHERE " + field_objectid + " = " 
        + lexical_cast<string>(object_id));

    Statement *stmt =
    _db.Prepare("INSERT INTO " + table + "(" + field_objectid + ", " + field_key + ", " + field_value
            + ") VALUES (?, ?, ?)");
    for(typename std::map<string, T>::iterator i = properties.begin();
            i!=properties.end(); ++i) {
        stmt->Bind(1, object_id);
        stmt->Bind(2, i->first);
        stmt->Bind(3, i->second);
        stmt->Step();
        stmt->Reset();
    }

    delete stmt;
}

template<typename T>
void StateSaverSQLite::ReadCustomProperties(int object_id, std::map<string, T> &properties,
        string table, const string field_objectid, const string field_key, const string field_value)
{
    Statement *stmt =
    _db.Prepare("SELECT " + field_key + ", " + field_value
            + " FROM " + table + " WHERE " + field_objectid + " = ?");

    stmt->Bind(1, object_id);
    while (stmt->Step() != SQLITE_DONE) {
        properties[stmt->Column<string>(0)] = stmt->Column<T>(1);
    }

    delete stmt;
}

