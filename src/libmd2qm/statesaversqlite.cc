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
    _db.Open(file.c_str(), SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE);
    _frame = 0;
    _qmtop=&qmtop;
}

void StateSaverSQLite::Close()
{
    _db.Close();
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

    stmt = _db.Prepare("INSERT INTO molecules (id, name) VALUES (?, ?)");

    int imol=0;
    for (MoleculeContainer::iterator iter = _qmtop->Molecules().begin();
            iter != _qmtop->Molecules().end() ; ++iter) {
        Molecule *mol=*iter;
        stmt->Bind(1, mol->getId());
        stmt->Bind(2, mol->getName());
        stmt->Step();
        stmt->Reset();
    }

    delete stmt;
}

void StateSaverSQLite::WriteCrgUnits(int frameid) {
    Statement *stmt = _db.Prepare(
            "INSERT INTO crgunits (id, name, energy, occ) VALUES (?,?,?,?)"
    );

    int imol=0;
    for (vector < QMCrgUnit *>::iterator iter = _qmtop->CrgUnits().begin(); iter!=_qmtop->CrgUnits().end(); ++iter) {
        QMCrgUnit *crg = *iter;
        stmt->Bind<int>(1, crg->getId());
        stmt->Bind(2, crg->getName());;
        stmt->Bind(3, crg->getEnergy());
        stmt->Bind(4, crg->getOccupationProbability());
        stmt->Step();
        stmt->Reset();
    }
    delete stmt;
}

void StateSaverSQLite::WriteBeads(int frameid) {
    Statement *stmt= _db.Prepare(
            "INSERT INTO beads (id,name,symmetry,type,resnr,mass,charge,crgunit,"
            "crgunit_index,pos_x,pos_y,pos_z,u_x,u_y,u_z,v_x,v_y,v_z,molid) "
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
        stmt->Bind<int>(8, bi->GetCrgUnit()->getId());;
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
            "INSERT INTO pairs (crgunit1, crgunit2, rate12, rate21, r_x, r_y,r_z)"
            " VALUES (?,?,?,?,?,?,?)");

    QMNBList &nblist = _qmtop->nblist();
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        QMPair *pair = *iter;
        QMCrgUnit *crg1 = (*iter)->first;
        QMCrgUnit *crg2 = (*iter)->second;

        stmt->Bind(1, (int)crg1->getId());
        stmt->Bind(2, (int)crg2->getId());
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
