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

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <exception>

#include <votca/csg/atomtable.h>

namespace votca { namespace csg {

using namespace std;

AtomTable::AtomTable(void){

    for(int ind=0;ind<this->Atoms.size();ind++){   
        Atom atm;
        atm.mass         = MassNumber.at(ind);
        atm.atomicNumber = AtomicNumber.at(ind);
        atomMap.insert(std::make_pair(Atoms.at(ind),atm));
    }
}

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

    double mass = atomMap[sym].mass;
    if(mass==0.0){
        string err_msg = "Invalid atom symbol in AtomTable::getMass. There is"
                         " no atom with symbol "+sym;
        throw invalid_argument(err_msg);
    }
    return mass;
}

int AtomTable::getAtomicNumber(string sym){
    double num = atomMap[sym].atomicNumber;
    if(num==0.0){
        string err_msg = "Invalid atom number in AtomTable::getAtomicNumber. "
                         "There is no atom with symbol "+sym;
        throw invalid_argument(err_msg);
    }
    return num;
}
}}
