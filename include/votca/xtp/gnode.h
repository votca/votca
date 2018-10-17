/*
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_XTP_GNODE_H
#define	VOTCA_XTP_GNODE_H

#include <votca/tools/vec.h>
#include <votca/xtp/glink.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/qmpair.h>
#include <vector>
#include <votca/xtp/huffmantree.h>


using namespace std;


namespace votca { namespace xtp {


class GNode
{
    public:
        GNode():occupied(false),occupationtime(0.0),escape_rate(0.0),hasdecay(false){};
        
        ~GNode(){};

        int id;
        bool occupied;
        bool injectable;
        double occupationtime;
        double escape_rate;
        bool hasdecay;
        Eigen::Vector3d position;
        std::vector<GLink> events;
        // stuff for Coulomb interaction:
        double siteenergy;
        double reorg_intorig; // UnCnN
        double reorg_intdest; // UcNcC
        void AddEvent(GNode* seg2, double rate12,const Eigen::Vector3d& dr, double Jeff2, double reorg_out);
        const double &getEscapeRate(){return escape_rate;}
        void InitEscapeRate();
        void AddDecayEvent(double decayrate);
        void ReadfromSegment(const Segment& seg, int carriertype);
        void AddEventfromQmPair(const QMPair& pair,int carriertype,std::vector<GNode>& nodes);
        
 
        GLink* findHoppingDestination(double p)const;
        void MakeHuffTree();

    private:
        huffmanTree<GLink> hTree;
        void organizeProbabilities(int id, double add);
        void moveProbabilities(int id);

};





}}

#endif	// VOTCA_XTP_GNODE_H

