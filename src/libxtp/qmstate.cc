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

#include <votca/xtp/qmstate.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <regex>
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace votca {
    namespace xtp {

    std::string QMStateType::ToString() const{
      std::string identifier="";
      switch (_type) {
        case QMStateType::Singlet: identifier = "S";
          break;
        case QMStateType::Triplet: identifier = "T";
          break;
        case QMStateType::PQPstate: identifier = "PQP";
          break;
        case QMStateType::DQPstate: identifier = "DQP";
          break;
        case QMStateType::KSstate: identifier = "KS";
          break;
        case QMStateType::Gstate: identifier = "N";
          break;
        case QMStateType::Hole: identifier = "H";
          break;
        case QMStateType::Electron: identifier = "E";
          break;
      }
      return identifier;
    }
    
    std::string QMStateType::ToLongString() const{
      std::string identifier="";
      switch (_type) {
        case QMStateType::Singlet: identifier = "singlet";
          break;
        case QMStateType::Triplet: identifier = "triplet";
          break;
        case QMStateType::PQPstate: identifier = "perturbative-quasiparticle";
          break;
        case QMStateType::DQPstate: identifier = "diagonalised-quasiparticle";
          break;
        case QMStateType::KSstate: identifier = "Kohn-Sham-orbital";
          break;
        case QMStateType::Gstate: identifier = "groundstate";
          break;
          case QMStateType::Hole: identifier = "hole";
          break;
        case QMStateType::Electron: identifier = "electron";
          break;
      }
      return identifier;
    }
    
    void QMStateType::FromString(const std::string& statetypestring){
      std::string lower = boost::algorithm::to_lower_copy(statetypestring);
      boost::trim(lower);
      if(lower=="s" || lower=="singlet"){
        _type=QMStateType::Singlet;
      }else if(lower=="t" || lower=="triplet" ){
        _type=QMStateType::Triplet;
      }else if(lower=="pqp" || lower=="perturbative-quasiparticle"){
        _type=QMStateType::PQPstate;
      }else if(lower=="dqp" || lower=="diagonalised-quasiparticle" || lower=="qpdiag"){
        _type=QMStateType::DQPstate;
      }else if(lower=="ks" || lower=="kohn-sham-orbital"){
        _type=QMStateType::KSstate;
      }else if(lower=="n" || lower=="groundstate" || lower=="gs"){
        _type=QMStateType::Gstate;
       }else if(lower=="h" || lower=="hole" ){
        _type=QMStateType::Hole;
       }else if(lower=="e" || lower=="electron"){
        _type=QMStateType::Electron;
      }else{
        throw std::runtime_error("Statetype:"+statetypestring+" not recognized");
      }
    }
    
    std::string QMState::ToLongString()const{
      int index=_index;
      if(_type==QMStateType::Singlet || _type==QMStateType::Triplet){
        index++;
      }else if(_type==QMStateType::Gstate){
        return _type.ToLongString();
      }
      std::string result=_type.ToLongString()+(boost::format(" %i") % index ).str();
      if(_transition){
        result="Groundstate to "+result;
      }
      return result;
    }
    
    std::string QMState::ToString()const{
      int index=_index;
      if(_type==QMStateType::Singlet || _type==QMStateType::Triplet){
        index++;
      }else if(_type==QMStateType::Gstate){
        return _type.ToString();
      }
      std::string result=_type.ToString()+(boost::format("%i") % index ).str();
      if(_transition){
        result="N2"+result;
      }
      return result;
    }
    
    
    int QMState::DetermineIndex(const std::string& statestring){
     
      std::smatch search;
      std::regex reg("[0-9]+");
      
      bool found_integer=std::regex_search(statestring,search,reg);
      if(!found_integer){
        throw std::runtime_error("Found no index in string: "+statestring);
      }
      if(search.size()>1){
        throw std::runtime_error("Found more than 1 index in string: "+statestring);
      }
      
      int index=boost::lexical_cast<int>(search.str(0));
        if(_type.isExciton() || _type==QMStateType::Electron || _type==QMStateType::Hole){
        index--;
      }
      return index;
    }
    
    
    QMStateType QMState::DetermineType(const std::string& statestring){
         std::regex reg("[^0-9]+");
         std::smatch search;
         
       bool found_typestring=std::regex_search(statestring,search,reg);
      if(!found_typestring){
        throw std::runtime_error("Found no type in string: "+statestring);
      }
      if(search.size()>1){
        throw std::runtime_error("Found more than one type in string: "+statestring);
      }
        QMStateType type;
        type.FromString(search.str(0));
        
        return type;
    }

    void QMState::FromString(const std::string& statestring){
      std::string lower = boost::algorithm::to_lower_copy(statestring);
      boost::trim(lower);
      std::string rest;
      if (boost::starts_with(lower, "n2")){
        _transition=true;
        rest=lower.substr(2);
      }else if(boost::starts_with(lower, "groundstate to")){
        _transition=true;
        rest=lower.substr(14);
      }else{
        rest=lower;
        _transition=false;
      }
      boost::trim(rest);
      
      _type=DetermineType(rest);
      if(_type!=QMStateType::Singlet && _transition==true){
          throw std::runtime_error("Transition states only exist for singlets.");
      }
      if(_type!=QMStateType::Gstate){
        _index=DetermineIndex(rest);
     }else{
          _index=-1;
     }
    }
  

    }
}
