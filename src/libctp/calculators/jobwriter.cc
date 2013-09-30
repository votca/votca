#include "jobwriter.h"
#include <votca/ctp/job.h>
#include <fstream>
#include <boost/format.hpp>

using boost::format;
namespace votca { namespace ctp {
    
    
void JobWriter::Initialize(Property *options) {   
    
    // REGISTER FUNCTIONS
    _key_funct["xqmultipole:ct"] = &JobWriter::xqmultipole_ct;
    _key_funct["xqmultipole:chrg"] = &JobWriter::xqmultipole_chrg;
    _key_funct["xqmultipole:kmc"] = &JobWriter::xqmultipole_kmc;    
    _key_funct["edft"] = &JobWriter::edft;
    _key_funct["idft"] = &JobWriter::idft;
    
    // SPLIT KEYS
    string keys = options->get("options.jobwriter.keys").as<string>();    
    Tokenizer tok_keys(keys, " ,\t\n");
    tok_keys.ToVector(_keys);
    
    // VALIDATE KEYS
    vector<string> ::iterator vsit;
    for (vsit = _keys.begin(); vsit != _keys.end(); ++vsit) {
        map<string,WriteFunct>::iterator it = _key_funct.find(*vsit);
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


bool JobWriter::EvaluateFrame(Topology *top) {
    
    vector<string> ::iterator vsit;
    for (vsit = _keys.begin(); vsit != _keys.end(); ++vsit) {
        map<string,WriteFunct>::iterator it = _key_funct.find(*vsit);
        if (it != _key_funct.end()) {
            cout << endl << "... ... " << it->first << flush;
            WriteFunct write = it->second;
            ((*this).*write)(top);
        }
    }
    return 0;
}


void JobWriter::xqmultipole_chrg(Topology *top) {
    
    // SET UP FILE STREAM
    ofstream ofs;
    string jobFile = "jobs_chrg.xml";    
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    int jobCount = 0;    
    vector<Segment*>::iterator sit1;
    
    // DEFINE PAIR CHARGE STATES
    vector<string > states;
    vector<string> ::iterator vit;
    states.push_back("n");
    states.push_back("e");
    states.push_back("h");

    
    // CREATE JOBS FOR ALL SEGMENTS AND STATES
    cout << endl;    
    for (sit1 = top->Segments().begin(); sit1 < top->Segments().end(); ++sit1) {
        Segment *seg1 = *sit1;

        int id1 = seg1->getId();
        string name1 = seg1->getName();
        
        for (vit = states.begin(); vit != states.end(); ++vit) {
            int id = ++jobCount;
            string s1 = *vit;
            string tag = (format("%1$d:%3$s:%2$s") % id1 % s1 % name1).str();                
            string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps")
                % id1 % name1 % s1).str();
            string stat = "AVAILABLE";

            Job job(id, tag, input, stat);
            job.ToStream(ofs,"xml");

            cout << "\r... ... # = " << jobCount << flush;
        }
    }
    
    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
}

    
void JobWriter::xqmultipole_kmc(Topology *top) {

    double cutoff = _options->get("options.jobwriter.cutoff").as<double>();
    
    // SET UP FILE STREAM
    ofstream ofs;
    string jobFile = "jobs_kmc.xml";    
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    int jobCount = 0;
    
    vector<Segment*>::iterator sit1;
    vector<Segment*>::iterator sit2;
    
    // DEFINE PAIR CHARGE STATES
    map< string, vector<string> > state1_state2;
    map< string, vector<string> > ::iterator mit;
    vector<string> ::iterator vit;
    state1_state2["n"] = vector<string>(1,"n");
    state1_state2["h"] = vector<string>(1,"h");
    state1_state2["h"].push_back("e");
    state1_state2["e"] = vector<string>(1,"e");
    state1_state2["e"].push_back("h");
    
    // CREATE JOBS FOR ALL PAIRS AND STATES
    cout << endl;    
    for (sit1 = top->Segments().begin(); sit1 < top->Segments().end(); ++sit1) {
        Segment *seg1 = *sit1;
        
        for (sit2 = sit1+1; sit2 < top->Segments().end(); ++sit2) {
            Segment *seg2 = *sit2;
            
            double dR = abs(top->PbShortestConnect(seg1->getPos(),seg2->getPos()));
            if (dR > cutoff) continue;
            
            int id1 = seg1->getId();
            string name1 = seg1->getName();
            int id2 = seg2->getId();
            string name2 = seg2->getName();        

            for (mit = state1_state2.begin(); mit != state1_state2.end(); ++mit) {
                for (vit = mit->second.begin(); vit != mit->second.end(); ++vit) {
                    int id = ++jobCount;
                    string s1 = mit->first;
                    string s2 = *vit;
                    string tag = (format("%1$d%2$s:%3$d%4$s:%5$1.2fnm")
                        % id1 % s1 % id2 % s2 % dR).str();                
                    string input = (format("%1$d:%2$s:MP_FILES/%2$s_%3$s.mps "
                        "%4$d:%5$s:MP_FILES/%5$s_%6$s.mps")
                        % id1 % name1 % s1 % id2 % name2 % s2).str();
                    string stat = "AVAILABLE";

                    Job job(id, tag, input, stat);
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


void JobWriter::xqmultipole_ct(Topology *top) {

    // SET UP FILE STREAM
    ofstream ofs;
    string jobFile = "jobs_ct.xml";    
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << endl;
    
    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    
    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    }
    
    // DEFINE PAIR CHARGE STATES
    map< string, vector<string> > state1_state2;
    map< string, vector<string> > ::iterator mit;
    vector<string> ::iterator vit;
    state1_state2["n"] = vector<string>(1,"n");
    state1_state2["h"] = vector<string>(1,"h");
    state1_state2["h"].push_back("e");
    state1_state2["e"] = vector<string>(1,"e");
    state1_state2["e"].push_back("h");
    
    // CREATE JOBS FOR ALL PAIRS AND STATES
    cout << endl;
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        
        int id1 = (*pit)->Seg1()->getId();
        string name1 = (*pit)->Seg1()->getName();
        int id2 = (*pit)->Seg2()->getId();
        string name2 = (*pit)->Seg2()->getName();        
        
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
                
                Job job(id, tag, input, stat);
                job.ToStream(ofs,"xml");
                
                cout << "\r... ... # = " << jobCount << flush;
            }
        }        
    }
    
    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
}
    
    
void JobWriter::edft(Topology *top) {

    string jobFile = "edft.jobs";   
    
    ofstream ofs;
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
 
    ofs << "<jobs>" << endl;   

    /* this is only good when ALL molecules shall be written out 
    for (vector<Segment*>::iterator sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
        int id = (*sit)->getId();
        string tag = (*sit)->getId();
        string input = "";
        string stat = "AVAILABLE";
        Job job(id, tag, input, stat);
        job.ToStream(ofs,"xml");
    }
    */

    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    

    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    } 

    map< int,Segment* > segments;
    map< int,Segment* >::iterator sit;

    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        
        int id1 = (*pit)->Seg1()->getId();
        int id2 = (*pit)->Seg2()->getId();
	segments[id1] = (*pit)->Seg1();
        segments[id2] = (*pit)->Seg2();

    }

    for (sit = segments.begin(); sit != segments.end(); ++sit) {
    
        int id = ++jobCount;
        
        string tag = (format("%1$s") % sit->first).str();
        string input = sit->second->getName();
        string stat = "AVAILABLE";
        Job job(id, tag, input, stat);
        job.ToStream(ofs,"xml");
    }
     

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
}
    
void JobWriter::idft(Topology *top) {

    string jobFile = "idft.jobs";   
    ofstream ofs;

    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
 
    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    

    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    }    

    ofs << "<jobs>" << endl;    
    
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        
        int id1 = (*pit)->Seg1()->getId();
        string name1 = (*pit)->Seg1()->getName();
        int id2 = (*pit)->Seg2()->getId();
        string name2 = (*pit)->Seg2()->getName();   

        int id = ++jobCount;
        string tag = (format("%1$s:%2$s") % id1 % id2 ).str(); 
        string input = (format("%1$s:%2$s") % name1 % name2 ).str();
        string stat = "AVAILABLE";
        Job job(id, tag, input, stat);
        job.ToStream(ofs,"xml");
    }

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
}


}}
