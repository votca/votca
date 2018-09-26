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

#ifndef VOTCA_XTP_CHARGECARRIER_H
#define	VOTCA_XTP_CHARGECARRIER_H

#include <votca/tools/vec.h>
#include <votca/xtp/gnode.h>



namespace votca { namespace xtp {
    
  
       
        
            class Chargecarrier
            {
                public:
                    Chargecarrier(): lifetime(0.0),steps(0) 
                    {
                        dr_travelled=Eigen::Vector3d::Zero();
                        node=NULL;
                    }
                    ~Chargecarrier(){};
                    bool hasNode(){return (node!=NULL);}
                    void updateLifetime(double dt) { lifetime+=dt;}
                    void updateOccupationtime(double dt) { node->occupationtime+=dt;}
                    void updateSteps(unsigned t) { steps+=t;}
                    void resetCarrier() { lifetime=0;steps=0; dr_travelled=Eigen::Vector3d::Zero();}
                    const double& getLifetime(){return lifetime;}
                    const unsigned& getSteps(){return steps;}
                    const int& getCurrentNodeId(){return node->id;}
                    double getCurrentEnergy(){return node->siteenergy;}
                    Eigen::Vector3d getCurrentPosition(){return node->position;}
                    double getCurrentEscapeRate(){return node->escape_rate;}
                    GNode * getCurrentNode(){return node;}
                    void settoNote(GNode *newnode){node=newnode;
                        node->occupied=true;}

                    void jumpfromCurrentNodetoNode(GNode *newnode){
                        node->occupied=false;
                        settoNote(newnode);
                    }
                    int id;
                    
                    Eigen::Vector3d dr_travelled;
                    
                private:
                    GNode *node;
                    double lifetime;
                    unsigned steps;
            };
     



}}


#endif	// VOTCA_XTP_CHARGECARRIER_H
