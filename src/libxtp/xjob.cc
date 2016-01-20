#include <votca/xtp/xjob.h>
#include <boost/format.hpp>


using boost::format;

namespace votca { namespace xtp {
    
    
    
XJob::XJob(int id, string tag, vector<Segment*> &qmSegs, 
     vector<string> &qmSegMps, Topology *top) :

     _id(id), _tag(tag), _top(top), _start_from_cpt(false),
     _qmSegs(qmSegs), _qmSegMps(qmSegMps), _ptop(NULL), _clean_ptop(true) {

     // Sanity checks
     assert(qmSegs.size() == qmSegMps.size());

     // Calc. center
     this->CalcCenterPos();
}


XJob::XJob(PolarTop *ptop, bool start_from_cpt)
        : _id(-1), _tag("__notag__"), _top(NULL), _start_from_cpt(start_from_cpt),
          _ptop(ptop), _clean_ptop(false) {
    
    _center = _ptop->getCenter();
    _isSegInCenter.clear();    
    vector<PolarSeg*>::iterator sit;
    for (sit = _ptop->QM0().begin(); sit < _ptop->QM0().end(); ++sit) {
        _isSegInCenter[(*sit)->getId()] = true;
    }
}


XJob::~XJob() {
    if (_clean_ptop) delete _ptop;
}


void XJob::setPolarTop(PolarTop *ptop) { 
    _ptop = ptop;
    setShellSizes(ptop->QM0().size(), ptop->MM1().size(), ptop->MM2().size());
}


void XJob::CalcCenterPos() {

    vec center = vec(0,0,0);
    double weight = 0;
    
    // Ref. point for PB correction
    vec refPt = vec(0,0,0);
    if (_qmSegs.size() == 0) { ; }
    else { refPt = _qmSegs[0]->getPos(); }

    // Calc. center
    for (unsigned int i = 0; i < _qmSegs.size(); ++i) {
         
        Segment *seg = _qmSegs[i];
        vec pbc_com = refPt + _top->PbShortestConnect(refPt, seg->getPos());
        
        center += seg->Atoms().size() * pbc_com;
        weight += seg->Atoms().size();
        
        this->_isSegInCenter[seg->getId()] = true;
    }
    
    _center = center / weight;   
}
    
    
void XJob::WriteInfoLine(FILE *out) {

   fprintf(out, "%5d %-20s  ",
                _id, _tag.c_str() );
   fprintf(out, "E_TOT %+4.7f ",
                _E_Tot );
   fprintf(out, "EPP %+4.7f "
                "EPU %+4.7f "
                "EUU %+4.7f ",
                _EPP,
                _EPU,
                _EUU);
   fprintf(out, "EF0_0 %+4.7f "
                "EF0_1 %+4.7f "
                "EF0_2 %+4.7f "
                "EF1_1 %+4.7f "
                "EF1_2 %+4.7f "
                "EM_0_ %+4.7f "
                "EM_1_ %+4.7f "
                "EM_2_ %+4.7f ",
                _EF_PAIR_PAIR,
                _EF_PAIR_SPH1,
                _EF_PAIR_SPH2,
                _EF_SPH1_SPH1, 
                _EF_SPH1_SPH2,
                _EM_PAIR,
                _EM_SPH1,      
                _EM_SPH2);
   fprintf(out, "ITER %3d SPHERE %4d SHELL %4d "
                "CENTER %4.7f %4.7f %4.7f ",
                _iter, _qm0_size+_mm1_size, _mm2_size, 
                _center.getX(), _center.getY(), _center.getZ() );
   fprintf(out, "E_CUT0_CUT0 %+4.7f "
                "E_CUT0_CUT1 %+4.7f "
                "E_CUT0_CUT2 %+4.7f "
                "E_CUT1_CUT1 %+4.7f "
                "E_CUT1_CUT2 %+4.7f "
                "E_PERM %+4.7f "
                "E_INDU %+4.7f ",
                _E_Pair_Pair,
                _E_Pair_Sph1,
                _E_Pair_Sph2,
                _E_Sph1_Sph1,  
                _E_Sph1_Sph2,
                _E_PERM,
                _E_INDU );
   fprintf(out, "\n");

}


Property XJob::GenerateOutputProperty() {
    
    Property prop;
    Property &out = prop.add("output","");
    Property *next = NULL;
    
    next = &out.add("summary", "");
    next->add("type", _tag);
    next->add("xyz", (format("%1$+1.7f %2$+1.7f %3$+1.7f") 
        % _center.getX() % _center.getY() % _center.getZ()).str())
        .setAttribute("unit","nm");
    next->add("total_scf", (format("%1$+1.7f") 
        % (_EF_PAIR_PAIR+_EF_PAIR_SPH1+_EF_PAIR_SPH2+_EF_SPH1_SPH1
          +_EF_SPH1_SPH2+_EM_PAIR+_EM_SPH1+_EM_SPH2)).str())
        .setAttribute("unit","eV");
    next->add("total", (format("%1$+1.7f") 
        % _E_Tot).str())
        .setAttribute("unit","eV");
    next->add("estat", (format("%1$+1.7f") 
        % _EPP).str())
        .setAttribute("unit","eV");
    next->add("eindu", (format("%1$+1.7f") 
        % _EPU).str())
        .setAttribute("unit","eV");
    
    next = &out.add("terms_i", "");
    next->add("F-00-01-02", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _EF_PAIR_PAIR % _EF_PAIR_SPH1 % _EF_PAIR_SPH2).str());
    next->add("F-11-12---", (format("%1$+1.5e %2$+1.5e") % _EF_SPH1_SPH1 % _EF_SPH1_SPH2).str());
    next->add("M-00-11-22", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _EM_PAIR % _EM_SPH1 % _EM_SPH2).str());
    next->add("E-PP-PU-UU", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _EPP % _EPU % _EUU).str());
    
    next = &out.add("shells", "");
    next->add("QM0", (format("%1$d") % _qm0_size).str());
    next->add("MM1", (format("%1$d") % _mm1_size).str());
    next->add("MM2", (format("%1$d") % _mm2_size).str());
    
    next = &out.add("convg", "");
    next->add("iter", (format("%1$d") % _iter).str());
        
    return prop;
}


void XJob::setInfoLine(bool printMM, bool printQM) {

    // Job ID & tag
    string str0 = ( format("%1$5d %2$-10s ") 
        % _id % _tag ).str();
    
    // Shell sizes & iterations
    string str1 = ( format("|QM0| %1$3d |MM1| %2$3d |MM2| %3$3d IT %4$2d ") 
        % _qm0_size % _mm1_size % _mm2_size % _iter).str();
    
    // Center coordinates
    string str2 = ( format("XYZ %1$+4.7f %2$+4.7f %3$+4.7f ") 
        % _center.getX() % _center.getY() % _center.getZ()).str();
    
    // QM energies
    string str3 = ( format("QMMM %3$+4.7e QM %1$+4.7e SF %2$+4.7e ") 
        % _E_QM % _E_SF % _E_QMMM ).str();
    
    // MM energies
    string str4 = ( format("TT %4$+4.7f PP %1$+4.7f PU %2$+4.7f UU %3$+4.7f ") 
        % _EPP % _EPU % _EUU % _E_Tot ).str();
    string str5 = ( format("F00 %1$+4.7f F01 %2$+4.7f F02 %3$+4.7f F11 %4$+4.7f F12 %5$+4.7f ") 
        % _EF_PAIR_PAIR % _EF_PAIR_SPH1 % _EF_PAIR_SPH2 % _EF_SPH1_SPH1 % _EF_SPH1_SPH2).str();    
    string str6 = ( format("M0 %1$+4.7f M1 %2$+4.7f M2 %3$+4.7f ") 
        % _EM_PAIR % _EM_SPH1 % _EM_SPH2 ).str();
    
    // Assemble
    _infoLine = str0;    
    if (printQM)
        _infoLine += str3;
    if (printMM)
        _infoLine += str4 + str5 + str6;
    _infoLine += str1 + str2;
    return;
}
    

template<typename JobContainer, typename pJob>
JobContainer XJOBS_FROM_TABLE(const string &job_file, Topology *top) {
    
    
    assert(false); // Define specialization
    JobContainer jobcnt;
    return jobcnt;
}


template<>
vector<Segment*> XJOBS_FROM_TABLE< vector<Segment*>, Segment* >(const string &job_file, Topology *top) {
    
    return top->Segments();
}


template<>
QMNBList XJOBS_FROM_TABLE< QMNBList, QMPair* >(const string &job_file, Topology *top) {
    
    return top->NBList();
}


template<>
vector<XJob*> XJOBS_FROM_TABLE< vector<XJob*>, XJob* >(const string &job_file, Topology *top) {
    
    vector<XJob*> xjobs;
    
    std::string line;
    std::ifstream intt;
    intt.open(job_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            vector<Segment*> qmSegs;
            vector<string>   qmSegMps;
            
            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

// Sample line
// # JOB_ID TAG    SEG1_ID SEG1_NAME SEG1_MPS SEG2_ID SEG2_NAME SEG2_MPS
//   1      E_CT   182:C60:MP_FILES/c60.mps   392:DCV:MP_FILES/dcv.mps
//   1  AVAILABLE E_CT 182:C60:MP_FILES/c60.mps 392:DCV:MP_FILES/dcv.mps

            int jobUserId   = boost::lexical_cast<int>(split[0]);
            string tag      = split[1];
            
            for (unsigned int i = 2; i < split.size(); ++i) {
                
                string id_seg_mps = split[i];
                vector<string> split_id_seg_mps;
                Tokenizer toker(id_seg_mps, ":");
                toker.ToVector(split_id_seg_mps);
                
                int segId = boost::lexical_cast<int>(split_id_seg_mps[0]);
                string segName = split_id_seg_mps[1];
                string mpsFile = split_id_seg_mps[2];
                
                Segment *seg = top->getSegment(segId);
                if (seg->getName() != segName) {
                    printf("\n\n");
                    printf("INPUT ERROR in file '%s'\n", job_file.c_str());
                    printf("Seg %d is named '%s', not '%s'\n", 
                            segId, seg->getName().c_str(), segName.c_str());
                    throw std::runtime_error("Input does not match topology.");
                }
                
                qmSegs.push_back(seg);
                qmSegMps.push_back(mpsFile);               
            }
            
            xjobs.push_back(new XJob(xjobs.size()+1, tag, qmSegs, qmSegMps, top));
            xjobs.back()->setUserId(jobUserId);

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << job_file << endl;
           throw runtime_error("Please supply input file.");           }
    
    return xjobs;
}


}}
