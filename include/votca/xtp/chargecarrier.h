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
 * author: Kordt
 */

#ifndef __VOTCA_CHARGECARRIER_H
#define	__VOTCA_CHARGECARRIER_H

// #include <votca/xtp/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath> // needed for abs(double)

#include <votca/tools/vec.h>
#include <votca/tools/matrix.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include <votca/tools/random2.h>


#include <votca/xtp/gnode.h>
#include <votca/ctp/qmcalculator.h>
using namespace std;

namespace votca { namespace xtp {
    
  
       
        
            class Chargecarrier
            {
                public:
                    Chargecarrier(): lifetime(0.0),steps(0) 
                    {
                        dr_travelled=tools::vec(0.0,0.0,0.0);
                        node=NULL;
                    }
                    ~Chargecarrier(){};
                    bool hasNode(){return (node!=NULL);}
                    void updateLifetime(double dt) { lifetime+=dt;}
                    void updateOccupationtime(double dt) { node->occupationtime+=dt;}
                    void updateSteps(unsigned t) { steps+=t;}
                    void resetCarrier() { lifetime=0;steps=0; dr_travelled=tools::vec(0.0,0.0,0.0);}
                    const double& getLifetime(){return lifetime;}
                    const unsigned& getSteps(){return steps;}
                    const int& getCurrentNodeId(){return node->id;}
                    double getCurrentEnergy(){return node->siteenergy;}
                    tools::vec getCurrentPosition(){return node->position;}
                    double getCurrentEscapeRate(){return node->escape_rate;}
                    GNode * getCurrentNode(){return node;}
                    void settoNote(GNode *newnode){node=newnode;
                        node->occupied=true;}

                    void jumpfromCurrentNodetoNode(GNode *newnode){
                        node->occupied=false;
                        settoNote(newnode);
                    }
                    int id;
                    
                    tools::vec dr_travelled;
                    
                private:
                    GNode *node;
                    double lifetime;
                    unsigned steps;
            };
     



}}


#endif	/* __VOTCA_CHARGECARRIER_H */
