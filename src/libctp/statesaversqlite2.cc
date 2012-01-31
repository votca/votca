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

#include <votca/ctp/statesaversqlite2.h>
#include <votca/tools/statement.h>

namespace votca { namespace ctp {

void StateSaverSQLite2::Open(Topology& qmtop, const string &file) {

    _sqlfile = file;
    _db.OpenHelper(file.c_str());
    
    _qmtop = &qmtop;
    _frames.clear();
    _topIds.clear();

    _frame = 0;
    _current_frame = -1;
    _was_read = false;

    // Query available frames in database
    Statement *stmt = _db.Prepare("SELECT _id, id FROM frames;");
    while(stmt->Step() != SQLITE_DONE) {
        _frames.push_back(stmt->Column<int>(0));
        _topIds.push_back(stmt->Column<int>(1));
    }
    delete stmt;
    
    
}



void StateSaverSQLite2::WriteFrame() {

    bool hasAlready = this->HasTopology(_qmtop);

    if ( ! hasAlready ) {
        if (_qmtop->getDatabaseId() >= 0 ) {
            throw runtime_error ("How was this topology generated? ");
        }
        _qmtop->setDatabaseId(_frames.size());
        _topIds.push_back(_qmtop->getDatabaseId());
        cout << "Saving ";
    }
    else {
        cout << "Updating ";
    }

    cout << "MD+QM topology ID " << _qmtop->getDatabaseId()
         << " (step = " << _qmtop->getStep()
         << ", time = " << _qmtop->getTime()
         << ") to " << _sqlfile << endl;
    cout << "... ";

    _db.BeginTransaction();    

    this->WriteMeta(hasAlready);
    this->WriteMolecules(hasAlready);
    this->WriteSegments(hasAlready);
    this->WriteFragments(hasAlready);
    this->WriteAtoms(hasAlready);

    _db.EndTransaction();

    cout << ". " << endl;
}


void StateSaverSQLite2::WriteMeta(bool update) {
    
    Statement *stmt;
    if ( update ) {
        return; // Nothing to do here
    }
    else {
        stmt = _db.Prepare( "INSERT INTO frames ("
                            "id,    time,  step,  "
                            "box11, box12, box13, "
                            "box21, box22, box23, "
                            "box31, box32, box33 )"
                            "VALUES ("
                            "?,     ?,     ?,"
                            "?,     ?,     ?,"
                            "?,     ?,     ?,"
                            "?,     ?,     ?)");
    }

    stmt->Bind(1, _qmtop->getDatabaseId());
    stmt->Bind(2, _qmtop->getTime());
    stmt->Bind(3, _qmtop->getStep());

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            stmt->Bind(4+3*i+j, _qmtop->getBox().get(i,j));
        }
    }
    stmt->InsertStep();
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::WriteMolecules(bool update) {
    cout << "Molecules";
    Statement *stmt;

    if (!update) {
        stmt = _db.Prepare("INSERT INTO molecules ("
                            "frame, top, id,"
                            "name, type    )"
                            "VALUES ("
                            "?,     ?,      ?,"
                            "?,     ?)");
    }
    else {
        return; // nothing to do here
    }
         

    stmt->Bind(1, _qmtop->getDatabaseId());

    vector < Molecule* > ::iterator mit;
    for (mit = _qmtop->Molecules().begin();
            mit < _qmtop->Molecules().end();
            mit++) {

        Molecule *mol = *mit;

        stmt->Bind(2, mol->getTopology()->getDatabaseId());
        stmt->Bind(3, mol->getId());
        stmt->Bind(4, mol->getName());
        stmt->Bind(5, mol->getName());

        stmt->InsertStep();
        stmt->Reset();
    }

    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::WriteSegments(bool update) {
    cout << ", segments";
    Statement *stmt;

    if (!update) {
        stmt = _db.Prepare("INSERT INTO segments ("
                            "frame, top, id,"
                            "name, type, mol)"
                            "VALUES ("
                            "?,     ?,  ?,"
                            "?,     ?,  ?)");
    }
    else {
        stmt = _db.Prepare("UPDATE segments "
                           "SET occ = ? "
                           "WHERE top = ? AND id = ?");
    }
    
    if (! update) { stmt->Bind(1, _qmtop->getDatabaseId()); }
    else { stmt->Bind(2, _qmtop->getDatabaseId()); }

    vector < Segment* > ::iterator sit;
    for (sit = _qmtop->Segments().begin();
            sit < _qmtop->Segments().end();
            sit++) {
        Segment *seg = *sit;

        if (!update) {

            stmt->Bind(2, seg->getTopology()->getDatabaseId());
            stmt->Bind(3, seg->getId());
            stmt->Bind(4, seg->getName());
            stmt->Bind(5, seg->getName());
            stmt->Bind(6, seg->getMolecule()->getId());
        }

        else {
            stmt->Bind(1, seg->getOcc());
            stmt->Bind(3, seg->getId());
        }

        stmt->InsertStep();
        stmt->Reset();

    }

    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::WriteFragments(bool update) {
    cout << ", fragments";

    Statement *stmt;

    if (! update) {
        stmt = _db.Prepare("INSERT INTO fragments ("
                            "frame, top, id,"
                            "name, type, mol,"
                            "seg, posX, posY,"
                            "posZ, symmetry )"
                            "VALUES ("
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?    )");
    }
    else {
        return; // nothing to do here
    }

    stmt->Bind(1, _qmtop->getDatabaseId());

    vector < Fragment* > ::iterator fit;
    for (fit = _qmtop->Fragments().begin();
            fit < _qmtop->Fragments().end();
            fit++) {
        Fragment *frag = *fit;

        stmt->Bind(2, frag->getTopology()->getDatabaseId());
        stmt->Bind(3, frag->getId());
        stmt->Bind(4, frag->getName());
        stmt->Bind(5, frag->getName());
        stmt->Bind(6, frag->getMolecule()->getId());
        stmt->Bind(7, frag->getSegment()->getId());
        stmt->Bind(8, frag->getPos().getX());
        stmt->Bind(9, frag->getPos().getY());
        stmt->Bind(10,frag->getPos().getZ());
        stmt->Bind(11,frag->getSymmetry());

        stmt->InsertStep();
        stmt->Reset();
    }

    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::WriteAtoms(bool update) {

    cout << ", atoms";

    Statement *stmt;
    if (! update) {
        stmt = _db.Prepare("INSERT INTO atoms ("
                            "frame, top, id,"
                            "name, type, mol,"
                            "seg, frag,  resnr,"
                            "resname, posX, posY,"
                            "posZ, weight)"
                            "VALUES ("
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?    )");
    }
    else {
        return; // nothing to do here
    }

    stmt->Bind(1, _qmtop->getDatabaseId());

    vector < Atom* > ::iterator ait;
    for (ait = _qmtop->Atoms().begin();
            ait < _qmtop->Atoms().end();
            ait++) {
        Atom *atm = *ait;
        
        stmt->Bind(2, atm->getTopology()->getDatabaseId());
        stmt->Bind(3, atm->getId());
        stmt->Bind(4, atm->getName());
        stmt->Bind(5, atm->getName());
        stmt->Bind(6, atm->getMolecule()->getId());
        stmt->Bind(7, atm->getSegment()->getId());
        stmt->Bind(8, atm->getFragment()->getId());
        stmt->Bind(9, atm->getResnr());
        stmt->Bind(10, atm->getResname());
        stmt->Bind(11, atm->getPos().getX());
        stmt->Bind(12, atm->getPos().getY());
        stmt->Bind(13, atm->getPos().getZ());
        stmt->Bind(14, atm->getWeight());

        stmt->InsertStep();
        stmt->Reset();
    }
    delete stmt;
    stmt = NULL;
}




bool StateSaverSQLite2::NextFrame() {

    _current_frame++;

    if(_current_frame >= _frames.size()) {
        return false;
    }
    else {
        this->ReadFrame();
        _was_read=true;
        return true;
    }

}


void StateSaverSQLite2::ReadFrame() {

    int topId = _topIds[_current_frame];

    cout << "Import MD+QM Topology ID " << topId
         << " (i.e. frame " << _current_frame << ")"
         << " from " << _sqlfile << endl;
    cout << "...";

    _qmtop->CleanUp();
    _qmtop->setDatabaseId(topId);

    this->ReadMeta(topId);
    this->ReadMolecules(topId);
    this->ReadSegments(topId);
    this->ReadFragments(topId);
    this->ReadAtoms(topId);

    cout << ". " << endl;
}


void StateSaverSQLite2::ReadMeta(int topId) {

    Statement *stmt = _db.Prepare("SELECT "
                                  "time, step, "
                                  "box11, box12, box13, "
                                  "box21, box22, box23, "
                                  "box31, box32, box33 "
                                  "FROM frames WHERE "
                                  "id = ?;");
    stmt->Bind(1, topId);

    if (stmt->Step() == SQLITE_DONE) {
        // ReadFrame should not have been called in the first place:
        throw runtime_error("Database appears to be broken. Abort...");
    }

    _qmtop->setTime(stmt->Column<double>(0));
    _qmtop->setStep(stmt->Column<int>(1));
    matrix boxv;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            boxv.set(i, j, stmt->Column<double>(2+3*i+j));
        }
    }
    _qmtop->setBox(boxv);
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::ReadMolecules(int topId) {

    cout << " Molecules";

    Statement *stmt = _db.Prepare("SELECT name "
                                  "FROM molecules "
                                  "WHERE top = ?;");
    stmt->Bind(1, topId);
    while (stmt->Step() != SQLITE_DONE) {
        Molecule *mol = _qmtop->AddMolecule(stmt->Column<string>(0));
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::ReadSegments(int topId) {

    cout << ", segments";

    Statement *stmt = _db.Prepare("SELECT name, mol, occ "
                                  "FROM segments "
                                  "WHERE top = ?;");
    stmt->Bind(1, topId);

    while (stmt->Step() != SQLITE_DONE) {

        string  name = stmt->Column<string>(0);
        int     id   = stmt->Column<int>(1);
        double  occ  = stmt->Column<double>(2);

        Segment *seg = _qmtop->AddSegment(name);
        seg->setOcc(occ);
        seg->setMolecule(_qmtop->getMolecule(id));
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::ReadFragments(int topId) {

    cout << ", fragments";

    Statement *stmt = _db.Prepare("SELECT "
                                  "name, mol, seg, "
                                  "posX, posY, posZ, "
                                  "symmetry "
                                  "FROM fragments "
                                  "WHERE top = ?;");

    stmt->Bind(1, topId);

    while (stmt->Step() != SQLITE_DONE) {

        string  name    = stmt->Column<string>(0);
        int     molid   = stmt->Column<int>(1);
        int     segid   = stmt->Column<int>(2);
        double  posX    = stmt->Column<double>(3);
        double  posY    = stmt->Column<double>(4);
        double  posZ    = stmt->Column<double>(5);
        int     symm    = stmt->Column<int>(3);

        Fragment *frag = _qmtop->AddFragment(name);
        frag->setSegment(_qmtop->getSegment(segid));
        frag->setMolecule(_qmtop->getMolecule(molid));
        frag->setPos(vec(posX, posY, posZ));
        frag->setSymmetry(symm);
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite2::ReadAtoms(int topId) {

    cout << ", atoms";

    Statement *stmt = _db.Prepare("SELECT "
                                  "name, mol, seg, frag, "
                                  "resnr, resname, "
                                  "posX, posY, posZ, "
                                  "weight "
                                  "FROM atoms "
                                  "WHERE top = ?;");

    stmt->Bind(1, topId);

    while (stmt->Step() != SQLITE_DONE) {

        string  name = stmt->Column<string>(0);
        int     molid = stmt->Column<int>(1);
        int     segid = stmt->Column<int>(2);
        int     fragid = stmt->Column<int>(3);
        int     resnr = stmt->Column<int>(4);
        string  resname = stmt->Column<string>(5);
        double  posX = stmt->Column<double>(6);
        double  posY = stmt->Column<double>(7);
        double  posZ = stmt->Column<double>(8);
        double  weight = stmt->Column<double>(9);

        Atom *atm = _qmtop->AddAtom(name);
        atm->setWeight(weight);
        atm->setPos( vec(posX, posY, posZ) );
        
        atm->setFragment(_qmtop->getFragment(fragid));
        atm->setSegment(_qmtop->getSegment(segid));
        atm->setMolecule(_qmtop->getMolecule(molid));

        atm->setResnr(resnr);
        atm->setResname(resname);  
    }
    delete stmt;
    stmt = NULL;
}





bool StateSaverSQLite2::HasTopology(Topology *top) {

    // Determine from topology ID whether database already stores a
    // (previous) copy

    Statement *stmt = _db.Prepare("SELECT id FROM frames");

    while (stmt->Step() != SQLITE_DONE) {
        if ( stmt->Column<int>(0) == top->getDatabaseId() ) { return true; }
        else { ; }
    }

    return false;

}


int StateSaverSQLite2::FramesInDatabase() {
    cout << "Reading file " << this->_sqlfile << ": Found " << _frames.size()
         << " frames stored in database. \n";
    return _frames.size();
}

}}
