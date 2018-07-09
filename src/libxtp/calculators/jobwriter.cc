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
    _key_funct["mps.single"] = &JobWriter::mps_single;
    _key_funct["mps.ct"] = &JobWriter::mps_ct;
    _key_funct["mps.chrg"] = &JobWriter::mps_chrg;
    _key_funct["mps.kmc"] = &JobWriter::mps_kmc;
    _key_funct["mps.background"] = &JobWriter::mps_background;

    
    // SPLIT KEYS
    string key="options."+Identify();
    std::string keys = options->get(key+".keys").as<string>();    
    Tokenizer tok_keys(keys, " ,\t\n");
    tok_keys.ToVector(_keys);
    
    // VALIDATE KEYS
    std::vector<std::string> ::iterator vsit;
    for (vsit = _keys.begin(); vsit != _keys.end(); ++vsit) {
        std::map<std::string,WriteFunct>::iterator it = _key_funct.find(*vsit);
        if (it == _key_funct.end()) {
            cout << endl << "... ... ERROR: Bad key '" << *vsit << "'" << endl;
            cout << endl << "Allowed keys: " << flush;
            for (it = _key_funct.begin(); it != _key_funct.end(); ++it)                
                cout << endl << "... " << it->first << flush;
            cout << endl << endl;
            throw runtime_error("No match for input key.");
        }
    }
    
    _options = options;        
    return;
}


bool JobWriter::EvaluateFrame(xtp::Topology *top) {
    
    std::vector<string> ::iterator vsit;
    for (vsit = _keys.begin(); vsit != _keys.end(); ++vsit) {
        std::map<string,WriteFunct>::iterator it = _key_funct.find(*vsit);
        if (it != _key_funct.end()) {
            cout << endl << "... ... " << it->first << flush;
            WriteFunct write = it->second;
            ((*this).*write)(top);
        }
    }
    return 0;
}


void JobWriter::mps_chrg(xtp::Topology *top) {
    
    // SET UP FILE STREAM
    ofstream ofs;
    std::string jobFile = Identify()+".mps.chrg.xml";
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    int jobCount = 0;    
    std::vector<xtp::Segment*>::iterator sit1;
    
    // DEFINE PAIR CHARGE STATES
    std::vector<string> states;
    std::vector<string> ::iterator vit;

    std::string str_states = _options->ifExistsReturnElseReturnDefault<string>("options.xjobwriter.states","n e h");
    std::string seg_pattern = _options->ifExistsReturnElseReturnDefault<string>("options.xjobwriter.pattern","*");
  
    Tokenizer tok_states(str_states, " ,\t\n");
    tok_states.ToVector(states);
    
    
    std::vector<string> seg_patterns;
    Tokenizer tok_pattern(seg_pattern, " ,\t\n");
    tok_pattern.ToVector(seg_patterns);

    
    // CREATE JOBS FOR ALL SEGMENTS AND STATES
    cout << endl;    
    for (sit1 = top->Segments().begin(); sit1 < top->Segments().end(); ++sit1) {
        xtp::Segment *seg1 = *sit1;

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
        
        for (vit = states.begin(); vit != states.end(); ++vit) {
            int id = ++jobCount;
            std::string s1 = *vit;
            std::string tag = (format("%1$d:%3$s:%2$s") % id1 % s1 % name1).str();                
            std::string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps")
                % id1 % name1 % s1).str();
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

    
void JobWriter::mps_kmc(xtp::Topology *top) {

    double cutoff = _options->get("options.xjobwriter.kmc_cutoff").as<double>();
    
    // SET UP FILE STREAM
    ofstream ofs;
    std::string jobFile = Identify()+".mps.kmc.xml";    
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    int jobCount = 0;
    
    std::vector<xtp::Segment*>::iterator sit1;
    std::vector<xtp::Segment*>::iterator sit2;
    
    // DEFINE PAIR CHARGE STATES
    std::map< string, vector<string> > state1_state2;
    std::map< string, vector<string> > ::iterator mit;
    std::vector<string> ::iterator vit;
    state1_state2["n"] = vector<string>(1,"n");
    state1_state2["h"] = vector<string>(1,"h");
    state1_state2["h"].push_back("e");
    state1_state2["e"] = vector<string>(1,"e");
    state1_state2["e"].push_back("h");
    
    // CREATE JOBS FOR ALL PAIRS AND STATES
    cout << endl;    
    for (sit1 = top->Segments().begin(); sit1 < top->Segments().end(); ++sit1) {
        xtp::Segment *seg1 = *sit1;
        
        for (sit2 = sit1+1; sit2 < top->Segments().end(); ++sit2) {
            xtp::Segment *seg2 = *sit2;
            
            double dR = abs(top->PbShortestConnect(seg1->getPos(),seg2->getPos()));
            if (dR > cutoff) continue;
            
            int id1 = seg1->getId();
            std::string name1 = seg1->getName();
            int id2 = seg2->getId();
            std::string name2 = seg2->getName();        

            for (mit = state1_state2.begin(); mit != state1_state2.end(); ++mit) {
                for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
                    int id = ++jobCount;
                    std::string s1 = mit->first;
                    std::string s2 = *vit;
                    std::string tag = (format("%1$d%2$s:%3$d%4$s:%5$1.2fnm")
                        % id1 % s1 % id2 % s2 % dR).str();                
                    std::string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps "
                        "%4$d:%5$s:MP_FILES/%5$s_%6$s.mps")
                        % id1 % name1 % s1 % id2 % name2 % s2).str();
                    std::string stat = "AVAILABLE";

                    xtp::Job job(id, tag, input, stat);
                    job.ToStream(ofs,"xml");
                    
                    cout << "\r... ... # = " << jobCount << flush;
                }
            }
        }
    }
    
    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
}


void JobWriter::mps_ct(xtp::Topology *top) {

    // SET UP FILE STREAM
    ofstream ofs;
    std::string jobFile = Identify()+".mps.ct.xml";
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    xtp::QMNBList::iterator pit;
    xtp::QMNBList &nblist = top->NBList();    
    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    }
    
    
    std::string str_states = _options->ifExistsReturnElseReturnDefault<string>("options.xjobwriter.states","nn eh");
    std::string seg_pattern = _options->ifExistsReturnElseReturnDefault<string>("options.xjobwriter.pattern","*");
  std::vector<string> states;
   
    Tokenizer tok_states(str_states, " ,\t\n");
    tok_states.ToVector(states);
    std::vector<string> seg_patterns;
    Tokenizer tok_pattern(seg_pattern, " ,\t\n");
    tok_pattern.ToVector(seg_patterns);
    
    
    
    // DEFINE PAIR CHARGE STATES
    std::map< string, vector<string> > state1_state2;
    std::map< string, vector<string> > ::iterator mit;
    std::vector<string> ::iterator vit;
    
    for (auto & state : states){
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
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        
        int id1 = (*pit)->Seg1()->getId();
        string name1 = (*pit)->Seg1()->getName();
        int id2 = (*pit)->Seg2()->getId();
        string name2 = (*pit)->Seg2()->getName(); 
        
        
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
        
        for (mit = state1_state2.begin(); mit != state1_state2.end(); ++mit) {
            for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
                int id = ++jobCount;
                string s1 = mit->first;
                string s2 = *vit;
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
    vector< xtp::Segment* > ::iterator sit;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
        ofs << (format("%1$4d %2$15s %3$-30s %4$-30s %5$-30s\n")
                % (*sit)->getId() 
                % (*sit)->getName()
                % ("MP_FILES/"+(*sit)->getName()+"_n.mps")
                % ("MP_FILES/"+(*sit)->getName()+"_e.mps")
                % ("MP_FILES/"+(*sit)->getName()+"_h.mps"));
    }
    ofs.close();
    return;
}


void JobWriter::mps_single(xtp::Topology *top) {
    
    // SET UP FILE STREAM
    ofstream ofs;
    std::string jobFile = Identify()+".mps.single.xml";
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    int jobCount = 0;    
   
    // DEFINE PAIR CHARGE STATES
    std::vector<std::string > states;
    std::vector<std::string> ::iterator vit;
    std::string str_states = _options->get("options.xjobwriter.states").as<std::string>();
    Tokenizer tok_states(str_states, " ,\t\n");
    tok_states.ToVector(states);
    
    // CREATE JOBS FOR ALL SEGMENTS AND STATES
    unsigned int single_id = _options->get("options.xjobwriter.single_id").as<int>();
    bool proceed = true;
    if (single_id < 1 || single_id > top->Segments().size()) {
        cout << endl 
             << "... ... ERROR Corrupt value in options.xjobwriter.single_id: "
             << "No such segment ID = " << single_id << ". Return." 
             << flush;
        ofs << "ERROR Corrupt value in options.xjobwriter.single_id" << endl;
        proceed = false;
    }
    
    if (proceed) {
        xtp::Segment *seg1 = top->getSegment(single_id);

        int id1 = seg1->getId();
        std::string name1 = seg1->getName();

        for (vit = states.begin(); vit != states.end(); ++vit) {
            int id = ++jobCount;
            std::string s1 = *vit;
            std::string tag = (format("%1$d:%3$s:%2$s") % id1 % s1 % name1).str();                
            std::string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps")
                % id1 % name1 % s1).str();
            std::string stat = "AVAILABLE";

            xtp::Job job(id, tag, input, stat);
            job.ToStream(ofs,"xml");
        }
    }

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    return;
}
    
    


}}
