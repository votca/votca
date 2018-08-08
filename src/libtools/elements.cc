/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/tools/elements.h>

using namespace votca::tools;
using namespace std;

/*************************
 * Public Facing Methods *
 *************************/

const double & Elements::getVdWChelpG(string name) const {
  if (_VdWChelpG.count(name) == 0)
    throw invalid_argument("Element not found in VdWChelpG map " + name);
  return _VdWChelpG.at(name);
}

const double & Elements::getVdWMK(string name) const {
  if (_VdWMK.count(name) == 0)
    throw invalid_argument("Element not found in VdWMP map " + name);
  return _VdWMK.at(name);
}

const double & Elements::getNucCrgECP(string name) const {
  throw invalid_argument("CrgECP map is deprecated");
}

const double & Elements::getPolarizability(string name) const {
  if (_ElPolarizability.count(name) == 0)
  throw invalid_argument("Element not found in ElPolarizability map " + name);
  return _ElPolarizability.at(name); 
}

double Elements::getCovRad(string name,string unit ) const { 
  //TODO - This should be replaced with an object, an object that should
  //       auto recognise the units and return it in a standard type
  if(!unit.compare("bohr")) return conv::ang2bohr*_CovRad.find(name)->second;
  if(!unit.compare("nm"))   return conv::ang2nm*_CovRad.find(name)->second;
  if(!unit.compare("ang"))  return _CovRad.find(name)->second; 
  throw invalid_argument("Must specify appropriate units " + 
      unit + " is not known");
}

const string & Elements::getEleName(int elenum) const {
  return _EleName.at(elenum);
}

const string & Elements::getEleShort(string elefull) const {
  return _EleShort.at(elefull);
}

bool  Elements::isMassAssociatedWithElement(double mass, double tolerance){
  auto closestMatch = findShortNameOfElementClosestInMass_(mass);
  if(closestMatch.second/_Mass[closestMatch.first]>tolerance) return false;
  return true;
}

string  Elements::getEleShortClosestInMass(double mass,double tolerance){
  auto closestMatch = findShortNameOfElementClosestInMass_(mass);
  if(closestMatch.second/_Mass[closestMatch.first]>tolerance){
    throw runtime_error("In attempt to determine if mass is associated "
        " with an element the mass exceeds tolerance of a possible match");
  }
  return closestMatch.first;
}

const string & Elements::getEleFull(string eleshort) const {
  return _EleFull.at(eleshort);
}

/*******************
 * Private Methods *
 *******************/

pair<string,double> Elements::findShortNameOfElementClosestInMass_(double mass){
  string eleShort = "H";
  double diff = abs(mass-_Mass[eleShort]);
  for(const auto & ele_pr : _Mass ){
    if(abs(ele_pr.second-mass)<diff){
      eleShort = ele_pr.first;
      diff = abs(ele_pr.second-mass);
    }
  }
  return pair<string,double>(eleShort,diff);
}

