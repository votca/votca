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

#include "jobwriter.h"
#include "votca/xtp/qmstate.h"
#include <votca/xtp/job.h>
#include <fstream>
#include <boost/format.hpp>
#include <votca/tools/tokenizer.h>

using boost::format;
namespace votca { namespace xtp {
    
    
void JobWriter::Initialize(Property *options) {   
    
    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );

    // REGISTER FUNCTIONS
    _key_funct["mps.dimer"] = &JobWriter::mps_dimer;
    _key_funct["mps.monomer"] = &JobWriter::mps_monomer;
    _key_funct["mps.background"] = &JobWriter::mps_background;
    
    // SPLIT KEYS
    string key="options."+Identify();
    std::string keys = options->get(key+".keys").as<string>();    
    Tokenizer tok_keys(keys, " ,\t\n");
    tok_keys.ToVector(_keys);
    
    // VALIDATE KEYS
    for (const std::string& key:_keys) {
        std::map<std::string,WriteFunct>::iterator it = _key_funct.find(key);
        if (it == _key_funct.end()) {
            cout << endl << "... ... ERROR: Bad key '" << key << "'" << endl;
            cout << endl << "Allowed keys: " << flush;
            for (auto& keyfunc:_key_funct)                
                cout << endl << "... " << keyfunc.first << flush;
            cout << endl << endl;
            throw runtime_error("No match for input key.");
        }
    }
    
    _options = options;        
    return;
}


bool JobWriter::EvaluateFrame(xtp::Topology *top) {
    
    for (const auto& key:_keys) {
            cout << endl << "... ... " << key<< flush;
            WriteFunct write = _key_funct.at(key);
            ((*this).*write)(top);
    }
    return 0;
}


void JobWriter::mps_monomer(xtp::Topology *top) {
    
    // SET UP FILE STREAM
    ofstream ofs;
    std::string jobFile = Identify()+".mps.monomer.xml";
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    int jobCount = 0;    

    std::string str_states = _options->ifExistsReturnElseReturnDefault<string>("options."+Identify()+".states","n e h");
    std::string seg_pattern = _options->ifExistsReturnElseReturnDefault<string>("options."+Identify()+".pattern","*");
  
    Tokenizer tok_states(str_states, " ,\t\n");
    std::vector<string> statesstring;
    tok_states.ToVector(statesstring);
    std::vector<QMState> states;
    for(const std::string& s:statesstring){
      states.push_back(QMState(s));
    }
    
    
    std::vector<string> seg_patterns;
    Tokenizer tok_pattern(seg_pattern, " ,\t\n");
    tok_pattern.ToVector(seg_patterns);

    
    // CREATE JOBS FOR ALL SEGMENTS AND STATES
    cout << endl;    
    for (xtp::Segment*seg1:top->Segments()) {
        int id1 = seg1->getId();
        std::string name1 = seg1->getName();

        bool do_continue=true;
        for (auto & pattern : seg_patterns){
            if (votca::tools::wildcmp(pattern.c_str(), name1.c_str())){
                do_continue=false;
                break;
            }
        }
        if(do_continue){continue;}
        
        for (QMState& state: states) {
            int id = ++jobCount;
            std::string tag = (format("%1$d:%3$s:%2$s") % id1 % state.ToString() % name1).str();                
            std::string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps")
                % id1 % name1 % state.ToString()).str();
            std::string stat = "AVAILABLE";
            xtp::Job job(id, tag, input, stat);
            job.ToStream(ofs,"xml");
            cout << "\r... ... # = " << jobCount << flush;
        }
    }
    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
}


void JobWriter::mps_dimer(xtp::Topology *top) {

    // SET UP FILE STREAM
    ofstream ofs;
    std::string jobFile = Identify()+".mps.dimer.xml";
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    xtp::QMNBList &nblist = top->NBList();    
    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    }
    
    
  std::string str_states = _options->ifExistsReturnElseReturnDefault<string>("options."+Identify()+".states","nn eh");
  std::string seg_pattern = _options->ifExistsReturnElseReturnDefault<string>("options."+Identify()+".pattern","*");
  std::vector<string> statesstring;
   
    Tokenizer tok_states(str_states, " ,\t\n");
    tok_states.ToVector(statesstring);
    std::vector<string> seg_patterns;
    Tokenizer tok_pattern(seg_pattern, " ,\t\n");
    tok_pattern.ToVector(seg_patterns);
    
    
    // DEFINE PAIR CHARGE STATES
    std::map< string, vector<string> > state1_state2;
    
    for (auto & state : statesstring){
        if(state=="nn"){
            state1_state2["n"] = vector<string>(1,"n");
        }
        else if ( state=="eh" || state=="he"){
            state1_state2["h"] = vector<string>(1,"e");
            state1_state2["e"] = vector<string>(1,"h");
        }
    }
    
    // CREATE JOBS FOR ALL PAIRS AND STATES
    cout << endl;
    for (xtp::QMPair* pair:nblist) {
        
        int id1 = pair->Seg1()->getId();
        string name1 = pair->Seg1()->getName();
        int id2 = pair->Seg2()->getId();
        string name2 = pair->Seg2()->getName(); 
        
        
        bool do_continue1=true;
        for (auto & pattern : seg_patterns){
            if (votca::tools::wildcmp(pattern.c_str(), name1.c_str())){
                do_continue1=false;
                break;
            }
        }
        bool do_continue2=true;
        for (auto & pattern : seg_patterns){
            if (votca::tools::wildcmp(pattern.c_str(), name2.c_str())){
                do_continue2=false;
                break;
            }
        }
        if(do_continue1 || do_continue2){continue;}
        
        for (auto& state1:state1_state2) {
            for (auto&  state2:state1.second) {
                int id = ++jobCount;
                string s1 = state1.first;
                string s2 = state2;
                string tag = (format("%1$d%2$s:%3$d%4$s") % id1 % s1 % id2 % s2).str();                
                string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps "
                        "%4$d:%5$s:MP_FILES/%5$s_%6$s.mps")
                        % id1 % name1 % s1 % id2 % name2 % s2).str();
                string stat = "AVAILABLE";
                
                xtp::Job job(id, tag, input, stat);
                job.ToStream(ofs,"xml");
                
                cout << "\r... ... # = " << jobCount << flush;
            }
        }        
    }
    
    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
}


void JobWriter::mps_background(xtp::Topology *top) {
    
    ofstream ofs;
    string tabFile = Identify()+".mps.background.tab";
    ofs.open(tabFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + tabFile);
    
    ofs << "# ID   TYPE    _n.mps    _e.mps    _h.mps \n";
    for (xtp::Segment* seg:top->Segments()) {
        ofs << (format("%1$4d %2$15s %3$-30s %4$-30s %5$-30s\n")
                % seg->getId() 
                % seg->getName()
                % ("MP_FILES/"+seg->getName()+"_n.mps")
                % ("MP_FILES/"+seg->getName()+"_e.mps")
                % ("MP_FILES/"+seg->getName()+"_h.mps"));
    }
    ofs.close();
    return;
}

}}
