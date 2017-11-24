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

#include <iostream>
#include <string>
#include <vector>

#include <votca/csg/atomtable.h>

namespace votca { namespace csg {

using namespace std;

int AtomTable::getIndex(string sym){
    int ind = 0;
    for(auto it=this->Atoms.begin();it!=this->Atoms.end();it++){
        if(*it==sym){
            return ind;
        }
        ind++;
    }
}

/* Check if a string is a legitimate atomic symbol */
bool AtomTable::checkSymbol(string sym){
  for(auto it=this->Atoms.begin();it!=this->Atoms.end();it++){
    if(*it==sym){
      return true;
    }
  }
  return false;
}

vector<string> AtomTable::getHalogens(void){
  vector<string> halogens = { "F","Cl","Br","I","At"};
  return halogens;
}

vector<string> AtomTable::getNoble(void){
  vector<string> noble = {"He","Ne","Ar","Kr","Xe","Rn","Uuo"};
  return noble;
}

double AtomTable::getMass(string sym){
    int ind = this->getIndex(sym);
    if(ind>=this->Atoms.size()){
        cerr << "Error Unrecognized symbol in atomtable.cc " << sym << endl;
        exit(1);
    } 

    double mass=this->MassNumber.at(ind);
    if(mass==-1.0){
        cerr << "Error mass for " << sym << " is unknown" << endl;
        exit(1);
    }
    return mass;
}

int AtomTable::getAtomicNumber(string sym){
    int ind = this->getIndex(sym);
    if(ind>=this->Atoms.size()){
        cerr << "Error Unrecognized symbol in atomtable.cc " << sym << endl;
        exit(1);
    } 

    return this->AtomicNumber.at(ind);
}
}}
