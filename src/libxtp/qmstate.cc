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

namespace votca {
    namespace xtp {

    std::string QMStateType::ToString() {
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
      }
      return identifier;
    }
    
    std::string QMStateType::ToLongString() {
      std::string identifier="";
      switch (_type) {
        case QMStateType::Singlet: identifier = "singlet state";
          break;
        case QMStateType::Triplet: identifier = "triplet state";
          break;
        case QMStateType::PQPstate: identifier = "perturbative quasiparticle state";
          break;
        case QMStateType::DQPstate: identifier = "diagonalised quasiparticle state";
          break;
        case QMStateType::KSstate: identifier = "Kohn Sham orbital";
          break;
        case QMStateType::Gstate: identifier = "Groundstate";
          break;
      }
      return identifier;
    }
    
    void QMStateType::FromString(const std::string& statetypestring){
      std::string lower = boost::algorithm::to_lower_copy(statetypestring);
      boost::trim(lower);
      if(lower=="s" || lower=="singlet" ){
        _type==QMStateType::Singlet;
      }else if(lower=="t" || lower=="triplet"){
        _type==QMStateType::Triplet;
      }else if(lower=="pqp" || lower=="perturbative quasiparticle state"){
        _type==QMStateType::PQPstate;
      }else if(lower=="dqp" || lower=="diagonalised quasiparticle state"){
        _type==QMStateType::DQPstate;
      }else if(lower=="ks" || lower=="kohn sham orbital"){
        _type==QMStateType::KSstate;
      }else if(lower=="n" || lower=="groundstate" || lower=="gs"){
        _type==QMStateType::Gstate;
      }else{
        throw runtime_error("Statetype:"+statetypestring+" not recognized");
      }
    }
    
    std::string QMState::ToLongString(){
      int index=_index;
      if(_type==QMStateType::Singlet || _type==QMStateType::Triplet){
        index++;
      }else if(_type==QMStateType::Gstate){
        return _type.ToLongString();
      }
      std::string result=_type.ToLongString()+(boost::format(" %i") % _index ).str();
      if(_transition){
        result="Groundstate to "+result;
      }
      return result;
    }
    
    std::string QMState::ToString(){
      int index=_index;
      if(_type==QMStateType::Singlet || _type==QMStateType::Triplet){
        index++;
      }else if(_type==QMStateType::Gstate){
        return _type.ToString();
      }
      std::string result=_type.ToString()+(boost::format("%i") % _index ).str();
      if(_transition){
        result="N2"+result;
      }
      return result;
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
      std::smatch integersearch;
      std::regex integer("(\\+|-)?[[:digit:]]+");
      
      bool found_integer=std::regex_search(rest,integersearch,integer);
      if(!found_integer){
        throw std::runtime_error("Found no index in string: "+rest);
      }
      if(integersearch.size()>1){
        throw std::runtime_error("Found more than 1 index in string: "+rest);
      }
      int index=0;
      for (auto& x:integersearch){
        index= boost::lexical_cast<int>(x);
      }
      
      
      
      
      if(_type==QMStateType::Singlet || _type==QMStateType::Triplet){
        index--;
      }
      _index=index;
    }
   
 
   
     
     
 
     
     
     
  

    }
}
