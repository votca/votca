/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_KMC_GNODE_H
#define	_VOTCA_KMC_GNODE_H

#include <votca/tools/vec.h>
#include <votca/xtp/glink.h>

using namespace std;
using namespace votca::xtp;


namespace votca { namespace xtp {

class GNode
{
    public:
        GNode():hasdecay(false){};
        ~GNode(){};

        int id;
        int occupied;
        int injectable;
        double occupationtime;
        double escaperate;
        bool hasdecay;
        votca::tools::vec position;
        vector<GLink> event;
        // stuff for Coulomb interaction:
        double siteenergy;
        double reorg_intorig; // UnCnN
        double reorg_intdest; // UcNcC
    
        void AddEvent(int seg2, double rate12, votca::tools::vec dr, double Jeff2, double reorg_out);
        const double &getEscapeRate(){return escaperate;}
        void InitEscapeRate();
        void AddDecayEvent(double _decayrate);
};





}}

#endif	/* _VOTCA_KMC_GNODE_H */

