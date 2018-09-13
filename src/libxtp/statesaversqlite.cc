/*
 *            Copyright 2009-2017 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/statesaversqlite.h>
#include <votca/tools/statement.h>

namespace votca { namespace xtp {

void StateSaverSQLite::Open(ctp::Topology& qmtop, const string &file, bool lock) {
    _sqlfile = file;
    if (lock) this->LockStateFile();
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
    if (lock) this->UnlockStateFile();
    return;    
}



void StateSaverSQLite::WriteFrame() {
    this->LockStateFile();
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
    this->WriteSegTypes(hasAlready);
    this->WriteSegments(hasAlready);
    this->WriteFragments(hasAlready);
    this->WriteAtoms(hasAlready);
    this->WritePairs(hasAlready);
    this->WriteSuperExchange(hasAlready);

    _db.EndTransaction();

    cout << ". " << endl;
    this->UnlockStateFile();
    return;
}


void StateSaverSQLite::WriteMeta(bool update) {
    
    Statement *stmt;
    if ( update ) {
        return; // Nothing to do here
    }
    else {
        stmt = _db.Prepare( "INSERT INTO frames ("
                            "id,    time,  step,  "
                            "box11, box12, box13, "
                            "box21, box22, box23, "
                            "box31, box32, box33, "
                            "canRigid )"
                            "VALUES ("
                            "?,     ?,     ?,"
                            "?,     ?,     ?,"
                            "?,     ?,     ?,"
                            "?,     ?,     ?,"
                            "?)");
    }

    stmt->Bind(1, _qmtop->getDatabaseId());
    stmt->Bind(2, _qmtop->getTime());
    stmt->Bind(3, _qmtop->getStep());

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            stmt->Bind(4+3*i+j, _qmtop->getBox().get(i,j));
        }
    }

    int canRigid = 0;
    if (_qmtop->canRigidify()) { canRigid = 1; }
    stmt->Bind(13, canRigid);

    stmt->InsertStep();
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::WriteMolecules(bool update) {
    cout << "Molecules" << flush;
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

    std::vector < ctp::Molecule* > ::iterator mit;
    for (mit = _qmtop->Molecules().begin();
            mit < _qmtop->Molecules().end();
            mit++) {

        ctp::Molecule *mol = *mit;

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


void StateSaverSQLite::WriteSegTypes(bool update) {
    cout << ", types" << flush;
    Statement *stmt;

    if (!update) {
        stmt = _db.Prepare("INSERT INTO segmentTypes ("
                           "frame, top, id,"
                           "name, basis, orbfile,"
                           "torbnrs, coordfile, canRigid)"
                           "VALUES ("
                           "?,  ?,  ?,"
                           "?,  ?,  ?,"
                           "?, ?,   ?)");
    }
    else {
        return; // nothing to do here
    }

    std::vector < ctp::SegmentType* > ::iterator typeit;
    for (typeit = _qmtop->SegmentTypes().begin();
            typeit < _qmtop->SegmentTypes().end();
            typeit++) {

        ctp::SegmentType *type = *typeit;

        if (!update) {
            stmt->Bind(1, _qmtop->getDatabaseId());
            stmt->Bind(2, type->getTopology()->getDatabaseId());
            stmt->Bind(3, type->getId());
            stmt->Bind(4, type->getName());
            stmt->Bind(5, type->getBasisName());
            stmt->Bind(6, type->getOrbitalsFile());            

            // Process transporting-orbital numbers
            // NOTE Obsolete, now read in only for IZindo calculator
            // string torbNrs = "";
            // std::vector <int> ::iterator iit;
            // for (int i = 0; i < type->getTOrbNrs().size(); i++) {
            //     int    nrInt = (type->getTOrbNrs())[i];
            //     string nrStr = boost::lexical_cast<string>(nrInt);
            //     torbNrs += ":" + nrStr;
            // }
            string torbStr = "NOT_USED";
            stmt->Bind(7, torbStr);
            stmt->Bind(8, type->getQMCoordsFile());

            int canRigid = 0;
            if (type->canRigidify()) { canRigid = 1; }
            stmt->Bind(9, canRigid);
        }
        stmt->InsertStep();
        stmt->Reset();
    }

    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::WriteSegments(bool update) {
    cout << ", segments" << flush;
    Statement *stmt;
    
     // Find out whether segments for this topology have already been created
    stmt = _db.Prepare("SELECT id FROM segments WHERE top = ?;");
    stmt->Bind(1, _qmtop->getDatabaseId());
    if (stmt->Step() == SQLITE_DONE) {        
        cout << " (create)" << flush;
    }
    else { 
        cout << " (recreate)" << flush;        
        Statement *stmt2 = _db.Prepare("DELETE FROM segments;");
        stmt2->Step();
        delete stmt2;
        Statement *stmt3 = _db.Prepare("UPDATE sqlite_sequence set seq = 0 where name='segments' ;");
        stmt3->Step();
        delete stmt3;
    }

    delete stmt;
    stmt = NULL;

    stmt = _db.Prepare("INSERT INTO segments ("
                        "frame, top, id,"
                        "name, type, mol,"
                        "posX, posY, posZ, "
                       "UnCnNe, UnCnNh, UcNcCe,"
                       "UcNcCh, UcCnNe, UcCnNh,"
                       "UnXnNs, UnXnNt, UxNxXs,"
                       "UxNxXt, UxXnNs, UxXnNt,"
                       "eAnion, eNeutral, eCation, eSinglet,eTriplet,"
                       "occPe, occPh,occPs, occPt,"
                       "has_e, has_h, has_s,has_t"
                       ") VALUES ("
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, ?, ?, ?, "
                           "?, ?, ?, ?, "
                           "?, ?, ?, ? "
                           ")");
   
    for (ctp::Segment *seg: _qmtop->Segments()) {
    
            stmt->Bind(1, _qmtop->getDatabaseId());
            stmt->Bind(2, seg->getTopology()->getDatabaseId());
            stmt->Bind(3, seg->getId());
            stmt->Bind(4, seg->getName());
            stmt->Bind(5, seg->getType()->getId());
            stmt->Bind(6, seg->getMolecule()->getId());
            stmt->Bind(7, seg->getPos().getX());
            stmt->Bind(8, seg->getPos().getY());
            stmt->Bind(9, seg->getPos().getZ());

            stmt->Bind(10, seg->getU_nC_nN(-1));
            stmt->Bind(11, seg->getU_nC_nN(+1));
            stmt->Bind(12, seg->getU_cN_cC(-1));
            stmt->Bind(13, seg->getU_cN_cC(+1));
            stmt->Bind(14, seg->getU_cC_nN(-1));
            stmt->Bind(15, seg->getU_cC_nN(+1));            
            stmt->Bind(16, seg->getU_nX_nN(+2));
            stmt->Bind(17, seg->getU_nX_nN(+3));
            stmt->Bind(18, seg->getU_xN_xX(+2));
            stmt->Bind(19, seg->getU_xN_xX(+3));          
            stmt->Bind(20, seg->getU_xX_nN(+2));
            stmt->Bind(21, seg->getU_xX_nN(+3));
            stmt->Bind(22, seg->getEMpoles(-1));
            stmt->Bind(23, seg->getEMpoles(0));
            stmt->Bind(24, seg->getEMpoles(1));
            stmt->Bind(25, seg->getEMpoles(2));
            stmt->Bind(26, seg->getEMpoles(3));
            stmt->Bind(27,seg->getOcc(-1));
            stmt->Bind(28,seg->getOcc(+1));
            stmt->Bind(29,seg->getOcc(+2));
            stmt->Bind(30,seg->getOcc(+3));
            
            int has_e = (seg->hasState(-1)) ? 1 : 0;
            int has_h = (seg->hasState(+1)) ? 1 : 0;
            int has_s = (seg->hasState(+2)) ? 1 : 0;
            int has_t = (seg->hasState(+3)) ? 1 : 0;
            stmt->Bind(31,has_e);
            stmt->Bind(32,has_h);
            stmt->Bind(33,has_s);
            stmt->Bind(34,has_t);
            
        stmt->InsertStep();
        stmt->Reset();

    }

    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::WriteFragments(bool update) {
    cout << ", fragments" << flush;

    Statement *stmt;

    if (! update) {
        stmt = _db.Prepare("INSERT INTO fragments ("
                            "frame, top, id,"
                            "name, type, mol,"
                            "seg, posX, posY,"
                            "posZ, symmetry, leg1,"
                            "leg2, leg3 )"
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

    for (ctp::Fragment *frag: _qmtop->Fragments()) {

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
        stmt->Bind(12,frag->getTrihedron()[0]);
        stmt->Bind(13,frag->getTrihedron()[1]);
        stmt->Bind(14,frag->getTrihedron()[2]);

        stmt->InsertStep();
        stmt->Reset();
    }

    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::WriteAtoms(bool update) {

    cout << ", atoms" << flush;

    Statement *stmt;
    if (! update) {
        stmt = _db.Prepare("INSERT INTO atoms ("
                            "frame, top, id,"
                            "name, type, mol,"
                            "seg, frag,  resnr,"
                            "resname, posX, posY,"
                            "posZ, weight, qmid,"
                            "qmPosX, qmPosY, qmPosZ,"
                            "element )"
                            "VALUES ("
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "?,     ?,  ?,"
                            "? )");
    }
    else {
        return; // nothing to do here
    }

    stmt->Bind(1, _qmtop->getDatabaseId());

    for (ctp::Atom *atm:_qmtop->Atoms()) {
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
        stmt->Bind(15, atm->getQMId());
        stmt->Bind(16, atm->getQMPos().getX());
        stmt->Bind(17, atm->getQMPos().getY());
        stmt->Bind(18, atm->getQMPos().getZ());
        stmt->Bind(19, atm->getElement());

        stmt->InsertStep();
        stmt->Reset();
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::WritePairs(bool update) {
    if ( ! _qmtop->NBList().size() ) { return; }
    
    cout << ", pairs" << flush;

    Statement *stmt;
    
    // Find out whether pairs for this topology have already been created
    stmt = _db.Prepare("SELECT id FROM pairs WHERE top = ?;");
    stmt->Bind(1, _qmtop->getDatabaseId());
    if (stmt->Step() == SQLITE_DONE) {        
        cout << " (create)" << flush;
    }
    else { 
        cout << " (recreate)" << flush;        
        Statement* stmt2 = _db.Prepare("DELETE FROM pairs;");
        stmt2->Step();
        delete stmt2;
        Statement*stmt3 = _db.Prepare("UPDATE sqlite_sequence set seq = 0 where name='pairs' ;");
        stmt3->Step();
        delete stmt3;
    }

    delete stmt;
    stmt = NULL;

     stmt = _db.Prepare("INSERT INTO pairs ("
                           "frame, top, id, "
                           "seg1, seg2, drX, "
                           "drY, drZ, "
                           "has_e, has_h,has_s,has_t, "
                           "lOe, lOh, lOs, lOt,"
                           "rate12e, rate21e, rate12h, rate21h,"
                           "rate12s, rate21s, rate12t, rate21t,"
                           "Jeff2e,  Jeff2h, Jeff2s, Jeff2t,"
                           "type "
                           ") VALUES ("
                           "?, ?, ?, "
                           "?, ?, ?, "
                           "?, ?, "
                           "?, ?, ?, ?, "
                           "?, ?, ?, ?, "
                           "?, ?, ?, ?, "
                           "?, ?, ?, ?, "
                           "?, ?, ?, ?, "
                           "? "
                           ")");

    for (ctp::QMPair *pair: _qmtop->NBList()) {

            int has_e = (pair->isPathCarrier(-1)) ? 1 : 0;
            int has_h = (pair->isPathCarrier(+1)) ? 1 : 0;
            int has_s = (pair->isPathCarrier(+2)) ? 1 : 0;
            int has_t = (pair->isPathCarrier(+3)) ? 1 : 0;

            stmt->Bind(1, _qmtop->getDatabaseId());
            stmt->Bind(2, pair->getTopology()->getDatabaseId());
            stmt->Bind(3, pair->getId());
            stmt->Bind(4, pair->Seg1PbCopy()->getId());
            stmt->Bind(5, pair->Seg2PbCopy()->getId());
            stmt->Bind(6, pair->R().getX());
            stmt->Bind(7, pair->R().getY());
            stmt->Bind(8, pair->R().getZ());
            stmt->Bind(9, has_e);
            stmt->Bind(10, has_h);
            stmt->Bind(11, has_s);
            stmt->Bind(12, has_t);
            stmt->Bind(13, pair->getLambdaO(-1));
            stmt->Bind(14, pair->getLambdaO(+1));
            stmt->Bind(15, pair->getLambdaO(+2));
            stmt->Bind(16, pair->getLambdaO(+3));
            stmt->Bind(17, pair->getRate12(-1));
            stmt->Bind(18, pair->getRate21(-1));
            stmt->Bind(19, pair->getRate12(+1));
            stmt->Bind(20, pair->getRate21(+1));
            stmt->Bind(21, pair->getRate12(+2));
            stmt->Bind(22, pair->getRate21(+2));
            stmt->Bind(23, pair->getRate12(+3));
            stmt->Bind(24, pair->getRate21(+3));
            stmt->Bind(25, pair->getJeff2(-1));
            stmt->Bind(26, pair->getJeff2(+1));
            stmt->Bind(27, pair->getJeff2(+2));
            stmt->Bind(28, pair->getJeff2(+3));
            stmt->Bind(29, (int)(pair->getType()) );
            
            stmt->InsertStep();
            stmt->Reset();
        }
    

    delete stmt;
    stmt = NULL;
}

void StateSaverSQLite::WriteSuperExchange(bool update) {
    if ( ! _qmtop->NBList().getSuperExchangeTypes().size() ) { return; }
    
    cout << ", super-exchange" << flush;

    Statement *stmt;
    
    // Find out whether pairs for this topology have already been created
    stmt = _db.Prepare("SELECT frame FROM superExchange WHERE top = ?;");
    stmt->Bind(1, _qmtop->getDatabaseId());
    if (stmt->Step() == SQLITE_DONE) {        
        cout << " (create)" << flush;
        delete stmt;
    }
    else { 
        cout << " (recreate)" << flush;        
        stmt = _db.Prepare("DELETE FROM superExchange;");
        stmt->Step();
        delete stmt;
        stmt = NULL;
        stmt = _db.Prepare("UPDATE sqlite_sequence set seq = 0 where name='pairs' ;");
        stmt->Step();
        delete stmt;
    }
    
    stmt = NULL;


    stmt = _db.Prepare("INSERT INTO superExchange ("
                           "frame, top, type"
                           ") VALUES ("
                           "?, ?, ?"
                           ")");

    list< ctp::QMNBList::SuperExchangeType* >::const_iterator seit;

    for (seit = _qmtop->NBList().getSuperExchangeTypes().begin();
         seit != _qmtop->NBList().getSuperExchangeTypes().end();
         seit++) {

        ctp::QMNBList::SuperExchangeType *seType = *seit;

        stmt->Bind(1, _qmtop->getDatabaseId());
        stmt->Bind(2, _qmtop->getDatabaseId());
        stmt->Bind(3, seType->asString());

        stmt->InsertStep();
        stmt->Reset();
    }

    delete stmt;
    stmt = NULL;    
    
    return;
}




bool StateSaverSQLite::NextFrame() {
    this->LockStateFile();
    bool hasNextFrame = false;
    _current_frame++;

    if(_current_frame < (int)_frames.size()) {
        this->ReadFrame();
        _was_read=true;
        hasNextFrame = true;
    }
    this->UnlockStateFile();
    return hasNextFrame;
}


void StateSaverSQLite::ReadFrame() {

    int topId = _topIds[_current_frame];

    cout << "Import MD+QM Topology ID " << topId
         << " (i.e. frame " << _current_frame << ")"
         << " from " << _sqlfile << endl;
    cout << "...";

    _qmtop->CleanUp();    
    _qmtop->setDatabaseId(topId);
    
    
    this->ReadMeta(topId);
    this->ReadMolecules(topId);
    this->ReadSegTypes(topId);
    this->ReadSegments(topId);
    this->ReadFragments(topId);
    this->ReadAtoms(topId);    
    this->ReadPairs(topId);
    this->ReadSuperExchange(topId);
    
    cout << ". " << endl;
}


void StateSaverSQLite::ReadMeta(int topId) {

    Statement *stmt = _db.Prepare("SELECT "
                                  "time, step, "
                                  "box11, box12, box13, "
                                  "box21, box22, box23, "
                                  "box31, box32, box33, "
                                  "canRigid "
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
    _qmtop->setCanRigidify(stmt->Column<int>(11));
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::ReadMolecules(int topId) {

    cout << " Molecules" << flush;

    Statement *stmt = _db.Prepare("SELECT name "
                                  "FROM molecules "
                                  "WHERE top = ?;");
    stmt->Bind(1, topId);
    while (stmt->Step() != SQLITE_DONE) {
	(void) _qmtop->AddMolecule(stmt->Column<string>(0));
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::ReadSegTypes(int topId) {

    cout << ", types";

    Statement *stmt = _db.Prepare("SELECT "
                                  "name, basis, orbfile, "
                                  "torbnrs, coordfile, canRigid "
                                  "FROM segmentTypes "
                                  "WHERE top = ?;");

    stmt->Bind(1, topId);
    while (stmt->Step() != SQLITE_DONE) {
        ctp::SegmentType *type = _qmtop->AddSegmentType(stmt->Column<string>(0));
        type->setBasisName(stmt->Column<string>(1));
        type->setOrbitalsFile(stmt->Column<string>(2));
        type->setQMCoordsFile(stmt->Column<string>(4));
        int canRigidify = stmt->Column<int>(5);
        type->setCanRigidify(canRigidify);

    }

    delete stmt;
    stmt = NULL;
}



void StateSaverSQLite::ReadSegments(int topId) {

    cout << ", segments" << flush;

    Statement *stmt = _db.Prepare("SELECT name, type, mol, "
                                  "posX, posY, posZ, "
                                  "UnCnNe, UnCnNh, UcNcCe,"
                                  "UcNcCh, UcCnNe , UcCnNh ,"            
                                  "UnXnNs, UnXnNt, UxNxXs ,"
                                  "UxNxXt , UxXnNs , UxXnNt,"
                                  "eAnion , eNeutral , eCation , eSinglet ,eTriplet ,"
                                  "occPe , occPh ,occPs , occPt ,"
                                  "has_e , has_h , has_s ,has_t "
                                  "FROM segments "
                                  "WHERE top = ?;");
    stmt->Bind(1, topId);

    while (stmt->Step() != SQLITE_DONE) {

        string  name = stmt->Column<string>(0);
        int     type = stmt->Column<int>(1);
        int     mId  = stmt->Column<int>(2);
        double  X    = stmt->Column<double>(3);
        double  Y    = stmt->Column<double>(4);
        double  Z    = stmt->Column<double>(5);
        
        double  l1   = stmt->Column<double>(6);
        double  l2   = stmt->Column<double>(7);
        double  l3   = stmt->Column<double>(8);
        double  l4   = stmt->Column<double>(9);
        double  l5   = stmt->Column<double>(10);
        double  l6   = stmt->Column<double>(11);
        
        double  x1   = stmt->Column<double>(12);
        double  x2   = stmt->Column<double>(13);
        double  x3   = stmt->Column<double>(14);
        double  x4   = stmt->Column<double>(15);
        double  x5   = stmt->Column<double>(16);
        double  x6   = stmt->Column<double>(17);
        
        double  e1   = stmt->Column<double>(18);
        double  e2   = stmt->Column<double>(19);       
        double  e3   = stmt->Column<double>(20);
        double  e4   = stmt->Column<double>(21);
        double  e5   = stmt->Column<double>(22);
        
        double  o1   = stmt->Column<double>(23);
        double  o2   = stmt->Column<double>(24);
        double  o3   = stmt->Column<double>(25);
        double  o4   = stmt->Column<double>(26);
        int     he   = stmt->Column<int>(27);
        int     hh   = stmt->Column<int>(28);
        int     hs   = stmt->Column<int>(29);
        int     ht   = stmt->Column<int>(30);

        bool has_e = (he == 1) ? true : false;
        bool has_h = (hh == 1) ? true : false;
        bool has_s = (hs == 1) ? true : false;
        bool has_t = (ht == 1) ? true : false;

        ctp::Segment *seg = _qmtop->AddSegment(name);
        seg->setMolecule(_qmtop->getMolecule(mId));
        seg->setType(_qmtop->getSegmentType(type));
        seg->setPos(vec(X, Y, Z));
        seg->setU_nC_nN(l1, -1);
        seg->setU_nC_nN(l2, +1);
        seg->setU_cN_cC(l3, -1);
        seg->setU_cN_cC(l4, +1);
        seg->setU_cC_nN(l5, -1);
        seg->setU_cC_nN(l6, +1);
        
        seg->setU_nX_nN(x1, +2);
        seg->setU_nX_nN(x2, +3);
        seg->setU_xN_xX(x3, +2);
        seg->setU_xN_xX(x4, +3);
        seg->setU_xX_nN(x5, +2);
        

        seg->setU_xX_nN(x6, +3);
                
        seg->setEMpoles(-1, e1);
        seg->setEMpoles(0, e2);
        seg->setEMpoles(1, e3);
        seg->setEMpoles(2, e4);
        seg->setEMpoles(3, e5);
        seg->setOcc(o1, -1);
        seg->setOcc(o2, +1);
        seg->setOcc(o3, +2);
        seg->setOcc(o4, +3);
        seg->setHasState(has_e, -1);
        seg->setHasState(has_h, +1);
        seg->setHasState(has_s, +2);
        seg->setHasState(has_t, +3);
        

        seg->getMolecule()->AddSegment(seg);
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::ReadFragments(int topId) {

    cout << ", fragments" << flush;

    Statement *stmt = _db.Prepare("SELECT "
                                  "name, mol, seg, "
                                  "posX, posY, posZ, "
                                  "symmetry, leg1, leg2, leg3 "
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
        int     symm    = stmt->Column<int>(6);
        int     leg1    = stmt->Column<int>(7);
        int     leg2    = stmt->Column<int>(8);
        int     leg3    = stmt->Column<int>(9);

        std::vector<int> trihedron;
        trihedron.push_back(leg1);
        if (leg2 >= 0) {trihedron.push_back(leg2);}
        if (leg3 >= 0) {trihedron.push_back(leg3);}

        ctp::Fragment *frag = _qmtop->AddFragment(name);
        frag->setSegment(_qmtop->getSegment(segid));
        frag->setMolecule(_qmtop->getMolecule(molid));
        frag->setPos(vec(posX, posY, posZ));
        frag->setSymmetry(symm);
        frag->setTrihedron(trihedron);

        frag->getSegment()->AddFragment(frag);
        frag->getMolecule()->AddFragment(frag);

    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::ReadAtoms(int topId) {

    cout << ", atoms" << flush;

    Statement *stmt = _db.Prepare("SELECT "
                                  "name, mol, seg, frag, "
                                  "resnr, resname, "
                                  "posX, posY, posZ, "
                                  "weight, qmid, qmPosX, "
                                  "qmPosY, qmPosZ, element "
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
        int     qmid = stmt->Column<int>(10);
        double  qmPosX = stmt->Column<double>(11);
        double  qmPosY = stmt->Column<double>(12);
        double  qmPosZ = stmt->Column<double>(13);
        string  element = stmt->Column<string>(14);

        ctp::Atom *atm = _qmtop->AddAtom(name);
        atm->setWeight(weight);
        atm->setQMPart(qmid, vec(qmPosX,qmPosY,qmPosZ));
        atm->setElement(element);
        atm->setPos( vec(posX, posY, posZ) );
        
        atm->setFragment(_qmtop->getFragment(fragid));
        atm->setSegment(_qmtop->getSegment(segid));
        atm->setMolecule(_qmtop->getMolecule(molid));

        atm->getFragment()->AddAtom(atm);
        atm->getSegment()->AddAtom(atm);
        atm->getMolecule()->AddAtom(atm);

        atm->setResnr(resnr);
        atm->setResname(resname);  
    }
    delete stmt;
    stmt = NULL;
}


void StateSaverSQLite::ReadPairs(int topId) {

    Statement *stmt = _db.Prepare("SELECT "
                                  "seg1, seg2,"
                                  "has_e, has_h,has_s,has_t,"
                                  "lOe, lOh, lOs, lOt,"
                                  "rate12e, rate21e, rate12h, rate21h,"
                                  "rate12s, rate21s, rate12t, rate21t,"
                                  "Jeff2e,  Jeff2h, Jeff2s, Jeff2t,"
                                  "type "
                                  "FROM pairs "
                                  "WHERE top = ?;");

    stmt->Bind(1, topId);
    ctp::QMNBList & nblist=_qmtop->NBList();

    while (stmt->Step() != SQLITE_DONE) {
        
        int     s1  = stmt->Column<int>(0);
        int     s2  = stmt->Column<int>(1);
        
        int     he  = stmt->Column<int>(2);
        int     hh  = stmt->Column<int>(3);
        int     hs  = stmt->Column<int>(4);
        int     ht  = stmt->Column<int>(5);
        double  l1  = stmt->Column<double>(6);
        double  l2  = stmt->Column<double>(7);
        double  l3  = stmt->Column<double>(8);
        double  l4  = stmt->Column<double>(9);
        double  r1  = stmt->Column<double>(10);
        double  r2  = stmt->Column<double>(11);
        double  r3  = stmt->Column<double>(12);
        double  r4  = stmt->Column<double>(13);
        double  r5  = stmt->Column<double>(14);
        double  r6  = stmt->Column<double>(15);
        double  r7  = stmt->Column<double>(16);
        double  r8  = stmt->Column<double>(17);
        double  je  = stmt->Column<double>(18);
        double  jh  = stmt->Column<double>(19);
        double  js  = stmt->Column<double>(20);
        double  jt  = stmt->Column<double>(21);
        int     tp  = stmt->Column<int>(22);
        
        ctp::QMPair *newPair = nblist.Add(_qmtop->getSegment(s1),
                                                _qmtop->getSegment(s2),false);
        
        bool has_e = (he == 0) ? false : true;
        bool has_h = (hh == 0) ? false : true;
        bool has_s = (hs == 0) ? false : true;
        bool has_t = (ht == 0) ? false : true;

        newPair->setIsPathCarrier(has_e, -1);
        newPair->setIsPathCarrier(has_h, +1);
        newPair->setIsPathCarrier(has_s, +2);
        newPair->setIsPathCarrier(has_t, +3);
        
        newPair->setLambdaO(l1, -1);
        newPair->setLambdaO(l2, +1);
        newPair->setLambdaO(l3, +2);
        newPair->setLambdaO(l4, +3);
        
        newPair->setRate12(r1, -1);
        newPair->setRate21(r2, -1);
        newPair->setRate12(r3, +1);
        newPair->setRate21(r4, +1);
        newPair->setRate12(r5, +2);
        newPair->setRate21(r6, +2);
        newPair->setRate12(r7, +3);
        newPair->setRate21(r8, +3);
        newPair->setJeff2(je, -1);
        newPair->setJeff2(jh, +1);
        newPair->setJeff2(js, +2);
        newPair->setJeff2(jt, +3);
        newPair->setType(tp);
    }
    cout << ", pairs" << flush;       
    delete stmt;
    stmt = NULL;

}


void StateSaverSQLite::ReadSuperExchange(int topId) {
    
    cout << ", super-exchange" << flush;

    Statement *stmt = _db.Prepare("SELECT "
                                  "type "
                                  "FROM superExchange "
                                  "WHERE top = ?;");

    stmt->Bind(1, topId);
    
    while (stmt->Step() != SQLITE_DONE) {
        string type = stmt->Column<string>(0);        
        _qmtop->NBList().AddSuperExchangeType(type);
    }
    delete stmt;
    stmt = NULL;    
    
    return;
}


bool StateSaverSQLite::HasTopology(ctp::Topology *top) {

    // Determine from topology ID whether database already stores a
    // (previous) copy

    Statement *stmt = _db.Prepare("SELECT id FROM frames");

    while (stmt->Step() != SQLITE_DONE) {
        if ( stmt->Column<int>(0) == top->getDatabaseId() ) { delete stmt;return true; }
        else { ; }
    }
    delete stmt;
    return false;

}


int StateSaverSQLite::FramesInDatabase() {
    cout << "Reading file " << this->_sqlfile << ": Found " << _frames.size()
         << " frames stored in database. \n";
    return _frames.size();
}


void StateSaverSQLite::LockStateFile() {    
    _flock = new boost::interprocess::file_lock(_sqlfile.c_str());
    _flock->lock();
    return;
}


void StateSaverSQLite::UnlockStateFile() {
    _flock->unlock();
    delete _flock;
    return;
}


}}
