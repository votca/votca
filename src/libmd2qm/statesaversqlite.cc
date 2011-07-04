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
    //ReadBeads();
//    ReadPairs();
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
        for(int j=0; j<9; ++j)
            m.set(i, j, stmt->Column<double>(2+i*3+j));
    _qmtop->setBox(m);
    delete stmt;
}

void StateSaverSQLite::ReadBeads() {
    Statement *stmt =
        _db.Prepare("SELECT "
            "rigidfrags._id, rigidfrags.id, rigidfrags.name, symmetry, type, resnr, mass, charge, molecule, "
            "conjseg_id, conjseg_index, pos_x, pos_y, pos_y, "
            "u_x, u_y, u_z, v_x, v_y, v_z "
            "FROM rigidfrags, molecules WHERE (molecules._id = molecule AND molecules.frame = ?)");
    stmt->Bind(1, _current_frame);

    while (stmt->Step() != SQLITE_DONE) {
/*
        byte_t symmetry =       ;
        string bead_name =      read<string> ();
        string type_name =      read<string> ();
        int resnr =             read<int>();
        // HACK: since only one residue is created this must be set to 0 by hand.
        resnr =0;
        double M =              read<double>();
        double Q =              read<double>();

        string crg_unit_name =  read<string> ();
        double energy = read<double>();

        unsigned short ipos =   read<unsigned short>();
        vec Pos =               read<vec> ();
        vec U =                 read<vec> ();
        vec V =                 read<vec> ();
        int molid =             read<int>();

        BeadType *type = _qmtop->GetOrCreateBeadType(type_name);
      
        QMBead *bead = dynamic_cast<QMBead*>(_qmtop->CreateBead(symmetry, bead_name, type, resnr, M, Q));
        _qmtop->getMolecule(molid)->AddBead(bead, bead_name);

        QMCrgUnit * acrg = _qmtop->GetCrgUnitByName(crg_unit_name);
        if(acrg == NULL)
            acrg = _qmtop->CreateCrgUnit(crg_unit_name, type_name, molid);
        acrg->setEnergy(energy);

        bead->setCrg(acrg);
        bead->setiPos(ipos);
        bead->setPos(Pos);
        bead->setU(U);
        bead->setV(V);
        bead->UpdateCrg();*/
    }
    delete stmt;
}

void StateSaverSQLite::WriteFrame()
{
    _db.BeginTransaction();
    Statement *stmt;
    stmt = _db.Prepare(
            "INSERT INTO frames (time, step, "
            "box11, box12, box13, box21, box22, box23, box31, box32, box33) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?)");

    stmt->Bind(1, _qmtop->getTime());
    stmt->Bind(2, _qmtop->getStep());

    for(int i=0; i<3; ++i)
        for(int j=0; j<9; ++j)
            stmt->Bind(3+i*3+j, _qmtop->getBox().get(i,j));

    stmt->Step();
    delete stmt;

    int frameid = _db.LastInsertRowId();
    WriteMolecules(frameid);
    WriteCrgUnits(frameid);
    WriteBeads(frameid);
    WritePairs(frameid);

    _db.EndTransaction();
    _frame++;
}

void StateSaverSQLite::WriteMolecules(int frameid)
{
    Statement *stmt;

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

void StateSaverSQLite::WriteCrgUnits(int frameid) {
    Statement *stmt = _db.Prepare(
            "INSERT INTO conjsegs (id, name, type, molecule, frame) VALUES (?,?,?,?,?)"
    );

    int imol=0;
    for (vector < QMCrgUnit *>::iterator iter = _qmtop->CrgUnits().begin(); iter!=_qmtop->CrgUnits().end(); ++iter) {
        QMCrgUnit *crg = *iter;
        stmt->Bind<int>(1, crg->getId());
        stmt->Bind(2, crg->getName());;
        stmt->Bind(3, crg->getType()->GetName());;
        stmt->Bind<int>(4, crg->getMolId());;
        stmt->Bind(5, frameid);;
        stmt->Step();
        stmt->Reset();
        crg->setDatabaseId(_db.LastInsertRowId());
    }
    delete stmt;
}

void StateSaverSQLite::ReadConjugatedSegments(void)
{
   Statement *stmt =
        _db.Prepare("SELECT name, type, molecule FROM conjsegs WHERE frame = ?;");
    stmt->Bind(1, _frames[_current_frame]);

    while(stmt->Step() != SQLITE_DONE) {
        QMCrgUnit *acrg = _qmtop->CreateCrgUnit(stmt->Column<string>(0), stmt->Column<string>(1), stmt->Column<int>(2));
    }
    delete stmt;
}

void StateSaverSQLite::WriteBeads(int frameid) {
    Statement *stmt= _db.Prepare(
            "INSERT INTO rigidfrags (id,name,symmetry,type,resnr,mass,charge,conjseg_id,"
            "conjseg_index,pos_x,pos_y,pos_z,u_x,u_y,u_z,v_x,v_y,v_z,molecule) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");

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
        stmt->Bind<int>(8, bi->GetCrgUnit()->getDatabaseId());;
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
        
        stmt->Step();
        stmt->Reset();
    }
}

void StateSaverSQLite::WritePairs(int frameid) {
    Statement *stmt;
    stmt = _db.Prepare(
            "INSERT INTO pairs (conjseg1, conjseg2, rate12, rate21, r_x, r_y,r_z)"
            " VALUES (?,?,?,?,?,?,?)");

    QMNBList &nblist = _qmtop->nblist();
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        QMPair *pair = *iter;
        QMCrgUnit *crg1 = (*iter)->first;
        QMCrgUnit *crg2 = (*iter)->second;

        stmt->Bind(1, (int)crg1->getDatabaseId());
        stmt->Bind(2, (int)crg2->getDatabaseId());
        stmt->Bind(3, pair->rate12());
        stmt->Bind(4, pair->rate21());
        stmt->Bind(5, pair->r().getX());
        stmt->Bind(6, pair->r().getY());
        stmt->Bind(7, pair->r().getZ());

        stmt->Step();
        stmt->Reset();
    }
    delete stmt;
}
