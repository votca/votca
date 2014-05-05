/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_GRAPHSQL_H_
#define __VOTCA_KMC_GRAPHSQL_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <votca/kmc/graph.h>
#include <votca/tools/database.h>

namespace votca { namespace kmc {

template<class TNode, class TLink>    
class GraphSQL : public Graph<TNode, TLink> {

public:
   
    /// Reads node information from filename
    void Initialize(string filename);
    
    bool &el_reorg() { return _el_reorg_fault; }    
    bool &ho_reorg() { return _ho_reorg_fault; }    

    bool _el_reorg_fault;
    bool _ho_reorg_fault;
};

template<class TNode, class TLink>    
inline void GraphSQL<TNode,TLink>::Initialize(string filename)
{   
    // Load Nodes
    votca::tools::Database db;
    db.Open( filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ, UnCnNe, UnCnNh, UcNcCe, UcNcCh, eAnion, eNeutral, eCation, UcCnNe, UcCnNh FROM segments;");

    while (stmt->Step() != SQLITE_DONE) 
    {    
        int id = stmt->Column<double>(0);
        double posX = stmt->Column<double>(1);
        double posY = stmt->Column<double>(2);
        double posZ = stmt->Column<double>(3);
        votca::tools::vec position(posX,posY,posZ);

        double UnCnNe = stmt->Column<double>(4);
        double UnCnNh = stmt->Column<double>(5);
        double UcNcCe = stmt->Column<double>(6);
        double UcNcCh = stmt->Column<double>(7);
        double eAnion = stmt->Column<double>(8);
        double eNeutral = stmt->Column<double>(9);
        double eCation = stmt->Column<double>(10);

        double UcCnNe = stmt->Column<double>(11);
        double UcCnNh = stmt->Column<double>(12);

        TNode* newTNode = this->AddNode(id,position);
        newTNode->setU(UnCnNe, UnCnNh, UcNcCe, UcNcCh);
        newTNode->setE(eAnion, eNeutral, eCation);
        newTNode->setu(UcCnNe, UcCnNh);
    }
    
    delete stmt;

    // Load Node Pairs

    stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', drX, drY, drZ, rate12e, rate12h, rate21e, rate21h, Jeff2e, Jeff2h, lOe, lOh  FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', -drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ', rate21e AS 'rate12e', rate21h AS 'rate12h', rate12e AS 'rate21e', rate12h AS 'rate21h',Jeff2e, Jeff2h, lOe, lOh  FROM pairs ORDER BY segment1;");
    _el_reorg_fault = false;
    _ho_reorg_fault = false;
    
    long id = 0;
    while (stmt->Step() != SQLITE_DONE) 
    {
        int node1_id = stmt->Column<int>(0);
        int node2_id = stmt->Column<int>(1);
        TNode* node1 = this->GetNode(node1_id);
        TNode* node2 = this->GetNode(node2_id);

        double drX = stmt->Column<double>(2);
        double drY = stmt->Column<double>(3);
        double drZ = stmt->Column<double>(4);
        votca::tools::vec r12(drX,drY,drZ);
        
        double rate12e = stmt->Column<double>(5);
        double rate12h = stmt->Column<double>(6);
        double rate21e = stmt->Column<double>(7);
        double rate21h = stmt->Column<double>(8);
        
        double Jeff2e = stmt->Column<double>(9);
        double Jeff2h = stmt->Column<double>(10);
        
        double lOe = stmt->Column<double>(11);
        double lOh = stmt->Column<double>(12);
        
        TLink* newTLink = this->AddLink(id,node1, node2, r12);
        newTLink->setRate(rate12e,rate12h,rate21e,rate21h);
        newTLink->setJeff2(Jeff2e,Jeff2h);
        newTLink->setlO(lOe,lOh);
        if(node1->UnCnNe() + node2->UcNcCe() + lOe == 0.0) { _el_reorg_fault = true;}
        if(node1->UnCnNh() + node2->UcNcCh() + lOh == 0.0) { _ho_reorg_fault = true;}
        //        newTLink->setcount(0);
        id++;
    }
        
    delete stmt;
    stmt = NULL;
}

}}



#endif

