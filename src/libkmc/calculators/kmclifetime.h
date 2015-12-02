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
 * author: Kordt
 */

#ifndef __VOTCA_KMC_LIFETIME_H
#define	__VOTCA_KMC_LIFETIME_H

// #include <votca/kmc/vssmgroup.h>
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
#include <votca/tools/statement.h>
#include <votca/tools/database.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include <votca/tools/random2.h>
#include <votca/kmc/kmccalculator.h>

#include <tr1/unordered_map>
#include <votca/kmc/gnode.h>
#include <math.h> // needed for fmod()

using namespace std;

namespace votca { namespace kmc {
    
   


class KMCLifetime : public KMCCalculator 
{
public:
    KMCLifetime() {};
   ~KMCLifetime() {};

    void Initialize(const char *filename, Property *options, const char *outputfile );
    bool EvaluateFrame();

private:
       
        
            class Chargecarrier
            {
                public:
                    Chargecarrier(): lifetime(0.0) 
                    {
                        dr_travelled=vec(0.0,0.0,0.0);
                    }
                    void updateLifetime(double dt) { lifetime+=dt;}
                    int position;
                    int id;
                    GNode *node;
                    vec dr_travelled;
                    
                private:
                    double lifetime;
            };
            
            
	    vector<GNode*>  LoadGraph();
            vector<double>  RunVSSM(vector<GNode*> node, unsigned int insertions, unsigned int numberofcharges, votca::tools::Random2 *RandomVariable);
            void InitialRates(vector<GNode*> node);
            void InitBoxSize(vector<GNode*> node);
            void progressbar(double fraction);
            void ReadLifetimeFile( string filename,vector<GNode*> node);
            double Promotetime(double cumulated_rate);
            void ResetForbidden(vector<int> &forbiddenid);
            void AddForbidden(int id, vector<int> &forbiddenid);
            bool Forbidden(int id, vector<int> forbiddenlist);
            //bool Surrounded(GNode* node, vector<int> forbiddendests);
            void printtime(int seconds_t);
            
            //vec _field;
            
            
            
            string _injection_name;
            string _injectionmethod;
           
            unsigned int _insertions;
            string _lifetimefile;
            double _maxrealtime;
            int _seed;
            int _numberofcharges;
            string _trajectoryfile;
            string _carriertype;
            double _temperature;
            string _outputfile;
            string _filename;
            string _rates;
};






}}


#endif	/* __VOTCA_KMC_MULTIPLE_H */