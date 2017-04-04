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
 * author: Kordt
 */

#ifndef __VOTCA_KMC_CALCULATOR_H
#define	__VOTCA_KMC_CALCULATOR_H

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
    
   


class KMCCalculator : public ctp::QMCalculator 
{
public:
    
    
   KMCCalculator();
   virtual ~KMCCalculator() {};
   
   
   
   virtual std::string  Identify() = 0;
   virtual void    Initialize(tools::Property *options) = 0;

   
protected:
       
        
            class Chargecarrier
            {
                public:
                    Chargecarrier(): lifetime(0.0),steps(0) 
                    {
                        dr_travelled=tools::vec(0.0,0.0,0.0);
                    }
                    ~Chargecarrier(){};
                    void updateLifetime(double dt) { lifetime+=dt;}
                    void updateSteps(unsigned t) { steps+=t;}
                    void resetCarrier() { lifetime=0;steps=0; dr_travelled=tools::vec(0.0,0.0,0.0);}
                    const double& getLifetime(){return lifetime;}
                    const unsigned& getSteps(){return steps;}
                    
                    int id;
                    GNode *node;
                    tools::vec dr_travelled;
                    
                private:
                    double lifetime;
                    unsigned steps;
            };
            
            std::string CarrierInttoLongString(int carriertype);
            std::string CarrierInttoShortString(int carriertype);
            int StringtoCarriertype(std::string name);
            
	    void LoadGraph(ctp::Topology *top);
            virtual void  RunVSSM(ctp::Topology *top){};
            void InitialRates();
            
            double Promotetime(double cumulated_rate);
            void ResetForbiddenlist(std::vector<int> &forbiddenid);
            void AddtoForbiddenlist(int id, std::vector<int> &forbiddenid);
            bool CheckForbidden(int id,const std::vector<int> &forbiddenlist);
            bool CheckSurrounded(GNode* node,const std::vector<int> &forbiddendests);
            GLink* ChooseHoppingDest(GNode* node);
            
            //tools::vec _field;
            
            std::vector<GNode*> _nodes;
            tools::Random2 * _RandomVariable;
           
            std::string _injection_name;
            std::string _injectionmethod;
            unsigned _outputsteps;
            
            int _seed;
            int _numberofcharges;
            tools::vec _field;
            int _carriertype;
            double _temperature;
            std::string _rates;
};






}}


#endif	/* __VOTCA_KMC_MULTIPLE_H */
