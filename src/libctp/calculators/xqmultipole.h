#ifndef XQMULTIPOLE_H
#define XQMULTIPOLE_H


#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/xinteractor.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>

namespace votca { namespace ctp {

class XQMP : public QMCalculator
{

public:

    XQMP() {};
   ~XQMP() {};

    string              Identify() { return "XQMultipole"; }
    void                Initialize(Topology *, Property *);
    bool                EvaluateFrame(Topology *top);

    // +++++++++++++++++++++++++++++++ //
    // MULTIPOLE ALLOCATION FROM INPUT //
    // +++++++++++++++++++++++++++++++ //

    void                Collect_JOB(string job_file, Topology *top);
    void                Collect_EMP(string emp_file, Topology *top);
    void                Collect_MPS(Topology *top);
    void                Collect_XML(string xml_file, Topology *top);

    void                Create_MPOLS(Topology *top);
    vector<APolarSite*>  Map_MPols_To_Seg(vector<APolarSite*> &, Segment *);
    vector<APolarSite*>  Parse_GDMA(string mps_file, int state);

    // +++++++++++++++++++++++++++ //
    // Job Operator (Thread class) //
    // +++++++++++++++++++++++++++ //

    class JobXQMP : public Thread
    {
    public:

        JobXQMP(int id,     Topology *top,        XQMP *master)
              : _id(id),         _top(top),     _master(master)
              {};

       ~JobXQMP() {};

        int         getId() { return _id; }
        void        setId(int id) { _id = id; }

        void        InitSlotData(Topology *top);
        void        Run(void);
        void        EvalJob(Topology *top, XJob *job);

    public:

        int                           _id;
        Topology                     *_top;
        XQMP                         *_master;
        XJob                         *_job;
        XInductor                     _inductor;

        vector< Segment* >            _segsPolSphere;  // Segments    in c/o 0-1
        vector< Segment* >            _segsOutSphere;  // Segments    in c/0 1-2
        vector< vector<APolarSite*> > _polsPolSphere;  // Polar sites in c/o 0-1
        vector< vector<APolarSite*> > _polsOutSphere;  // Polar sites in c/o 1-2
        vector< vector<APolarSite*> > _polarSites;     // Copy of top polar sites
        vector< vector<APolarSite*> > _polarSites_job; // Adapted to job specs

    };

    XJob               *RequestNextJob(int id, Topology *top);
    void                LockCout() { _coutMutex.Lock(); }
    void                UnlockCout() { _coutMutex.Unlock(); }


private:

    // ++++++++++++++++++++++++++++++++++++++++ //
    // MULTIPOLE ALLOCATION TO SEGMENTS / PAIRS //
    // ++++++++++++++++++++++++++++++++++++++++ //

    string                         _job_file;
    string                         _emp_file;
    string                         _xml_file;

    // Job info : JOB_ID PAIR_ID MPS_1 MPS_2 TAG
    vector<XJob*>                  _XJobs;
    vector<XJob*>::iterator        _nextXJob;
    Mutex                          _nextJobMutex;
    Mutex                          _coutMutex;
    Mutex                          _logMutex;
    bool                           _maverick;

    vector<int>                    _jobIds;
    map<int,string>                _jobId_jobTag;
    map<int,int>                   _jobId_pairId;
    map<int,string>                _jobId_mpsFile1;
    map<int,string>                _jobId_mpsFile2;
    map<int, pair<int,int> >       _pairId_seg1Id_seg2Id;
    vector<int>::iterator          _nextJob;

    // Allocate .mps files to segments (n 0, e -1, h +1)
    map<int,string>                 _segId_mpsFile_n;
    map<int,string>                 _segId_mpsFile_e;
    map<int,string>                 _segId_mpsFile_h;

    // Store polar site containers, one for each .mps file
    map<string,vector<APolarSite*> > _mpsFile_pSites;
    map<string,vector<APolarSite*> > _mpsFile_pSites_job;

    // Allocate polar sites to fragments in segments
    map<string, bool>               _map2md;
    map<string, vector<int> >       _alloc_frag_mpoleIdx;
    map<string, vector<string> >    _alloc_frag_mpoleName;
    map<string, vector<int> >       _alloc_frag_trihedron;
    map<string, vector<double> >    _alloc_frag_weights;

    string                          _pdb_check;
    bool                            _write_chk;
    string                          _write_chk_suffix;
    bool                            _chk_split_dpl;
    double                          _chk_dpl_spacing;
    string                          _chk_format;


    // ++++++++++++++++++++++++++++++++++++++++ //
    // INDUCTION + ENERGY EVALUATION            //
    // ++++++++++++++++++++++++++++++++++++++++ //

    // Control options
    bool                            _induce;
    bool                            _induce_intra_pair;
    int                             _subthreads;

    // Interaction parameters
    bool                            _useCutoff;
    double                          _cutoff;
    double                          _cutoff2;
    bool                            _useExp;
    double                          _aDamp;
    bool                            _useScaling;
    vector<double>                  _scale1;

    // Convergence parameters
    float                           _wSOR_N;
    float                           _wSOR_C;
    double                          _epsTol;
    int                             _maxIter;

    // Logging
    string                          _outFile;
    bool                            _energies2File;
    map<int,vector<int> >           _log_seg_iter;
    map<int,int>                    _log_seg_sphereSize;
    map<int,vec>                    _log_seg_com;


};

// ========================================================================== //
//                        XMULTIPOLE MEMBER FUNCTIONS                         //
// ========================================================================== //


void XQMP::Initialize(Topology *top, Property *opt) {

    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;
    

    string key = "options.xqmultipole.multipoles";

        if ( opt->exists(key) ) {
            _xml_file = opt->get(key).as< string >();
        }

    key = "options.xqmultipole.control";

        if ( opt->exists(key+".job_file")) {
            _job_file   = opt->get(key+".job_file").as<string>();
        }
        else {
            _job_file   = opt->get(key+".job_file").as<string>();
        }

        if ( opt->exists(key+".emp_file")) {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }
        else {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }

        if ( opt->exists(key+".output") ) {
            _outFile = opt->get(key+".output").as< string >();
            _energies2File = true;
        }
        else { _energies2File = false; }

        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<string>();
        }
        else { _pdb_check = ""; }

        if (opt->exists(key+".write_chk")) {
            _write_chk_suffix = opt->get(key+".write_chk").as<string>();
            _write_chk = true;
        }
        else { _write_chk = false; }

        if (opt->exists(key+".format_chk")) {
            _chk_format = opt->get(key+".format_chk").as<string>();
        }
        else { _chk_format = "xyz"; }

        if (opt->exists(key+".split_dpl")) {
            _chk_split_dpl = (opt->get(key+".split_dpl").as<int>() == 1) ?
                         true : false;
        }
        else { _chk_split_dpl = true; }

        if (opt->exists(key+".dpl_spacing")) {
            _chk_dpl_spacing = opt->get(key+".dpl_spacing").as<double>();
        }
        else {
            _chk_dpl_spacing = 1.0e-6;
        }


    key = "options.xqmultipole.tholemodel";

        if ( opt->exists(key+".induce") ) {
            int induce = opt->get(key+".induce").as< int >();
            _induce = (induce == 0) ? false : true;
        }
        else { _induce = true; }

        if ( opt->exists(key+".induce_intra_pair") ) {
            int induce = opt->get(key+".induce_intra_pair").as< int >();
            _induce_intra_pair = (induce == 0) ? false : true;
        }
        else { _induce_intra_pair = true; }

        if ( opt->exists(key+".cutoff1") ) {
            _cutoff = opt->get(key+".cutoff1").as< double >();
            if (_cutoff) { _useCutoff = true; }
        }
        if ( opt->exists(key+".cutoff2") ) {
            _cutoff2 = opt->get(key+".cutoff2").as< double >();
        }
        else {
            _cutoff2 = _cutoff;
        }
        if ( opt->exists(key+".subthreads") ) {
            _subthreads = opt->get(key+".subthreads").as< double >();
        }
        else {
            _subthreads = 1;
        }
        if ( opt->exists(key+".exp_damp") ) {
            _aDamp = opt->get(key+".exp_damp").as< double >();
            if (_aDamp) { _useExp = true; }
        }
        else {
            cout << endl << "... ... WARNING: No sharpness parameter supplied";
            cout << endl << "... ... ... Using default a = 0.39";
            _aDamp = 0.39;
            _useExp = true;
        }
        if ( opt->exists(key+".scaling") ) {
            _scale1 = opt->get(key+".scaling").as< vector<double> >();
            if (0 < _scale1.size() && _scale1.size() < 4) {
                _useScaling = true; }
            else {
                _useScaling = false;
                cout << endl << "... ... WARNING: 1-N SCALING SWITCHED OFF"; }
        }

    key = "options.xqmultipole.convergence";

        if ( opt->exists(key+".wSOR_N") ) {
            _wSOR_N = opt->get(key+".wSOR_N").as< float >();
        }
        else { _wSOR_N = 0.75; }
        if ( opt->exists(key+".wSOR_C") ) {
            _wSOR_C = opt->get(key+".wSOR_C").as< float >();
        }
        else { _wSOR_C = 0.75; }

        if ( opt->exists(key+".max_iter") ) {
            _maxIter = opt->get(key+".max_iter").as< int >();
        }
        else { _maxIter = 512; }

        if ( opt->exists(key+".tolerance") ) {
            _epsTol = opt->get(key+".tolerance").as< double >();
        }
        else { _epsTol = 0.001; }    

}


void XQMP::Collect_XML(string xml_file, Topology *top) {

    cout << endl 
         << "... ... ... Allocate polar sites to fragments. "
         << flush;

    string allocFile = xml_file;

    // ++++++++++++++++++++++++++++++++ //
    // Load polar-site indices from XML //
    // ++++++++++++++++++++++++++++++++ //

    // => Output to maps:
    map<string, vector<int> > alloc_frag_mpoleIdx;
    map<string, vector<string> > alloc_frag_mpoleName;
    map<string, vector<int> > alloc_frag_trihedron;
    map<string, vector<double> > alloc_frag_weights;

    Property allocation; // <- Which polar sites are part of which fragment?
    load_property_from_xml(allocation, allocFile.c_str());


    /* --- MULTIPOLES.XML Structure ---
     *
     * <topology>
     *
     *     <molecules>
     *          <molecule>
     *          <name></name>
     *
     *          <segments>
     *
     *              <segment>
     *              <name>DCV</name>
     *
     *              <map2md></map2md>
     *
     *              <fragments>
     *                  <fragment>
     *                  <name></name>
     *                  <mpoles></mpoles>
     *                  </fragment>
     *              </fragments>
     *              ...
     *              ...
     */


    string key = "topology.molecules.molecule";
    list<Property *> mols = allocation.Select(key);
    list<Property *>::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); molit++) {

        string molName = (*molit)->get("name").as<string> ();

        key = "segments.segment";
        list<Property *> segs = (*molit)->Select(key);
        list<Property *>::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); segit++) {

            string segName = (*segit)->get("name").as<string> ();

            // Default: Project multipoles onto rigidified coordinates
            if ( (*segit)->exists("map2md")) {
                int map2md = (*segit)->get("map2md").as<int>();
                _map2md[segName] = (map2md == 0) ? false : true;
            }
            else {
                _map2md[segName] = false; // i.e. map to rigidified coordinates
            }

            key = "fragments.fragment";
            list<Property *> frags = (*segit)->Select(key);
            list<Property *>::iterator fragit;

            for (fragit = frags.begin(); fragit != frags.end(); fragit++) {

                string fragName = (*fragit)->get("name").as<string> ();
                string mapKeyName = fragName + segName + molName;

                string mpoles = (*fragit)->get("mpoles").as<string> ();

                // Local frame for polar sites
                vector<int> trihedron_mps;
                if ((*fragit)->exists("localframe_mps")) {
                   cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Defining distinct local frame for polar sites."
                         << flush;
                   trihedron_mps = (*fragit)->get("localframe_mps")
                                         .as< vector<int> >();
                }
                else {
                   trihedron_mps = (*fragit)->get("localframe")
                                         .as< vector<int> >();
                }
                
                // Mapping weights for polar sites
                vector<double> weights_mps;
                if ((*fragit)->exists("weights_mps")) {
                    cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Using distinct weights for polar sites."
                         << flush;
                   weights_mps = (*fragit)->get("weights_mps")
                                       .as< vector<double> >();
                }
                else {
                   weights_mps = (*fragit)->get("weights")
                                       .as< vector<double> >();
                }

                Tokenizer tokPoles(mpoles, " \t\n");
                vector<string> mpoleInfo;
                tokPoles.ToVector(mpoleInfo);

                vector<int> mpoleIdcs;
                vector<string> mpoleNames;

                vector<string> ::iterator strit;
                for (strit=mpoleInfo.begin(); strit<mpoleInfo.end(); strit++) {

                    Tokenizer tokPoleInfo( (*strit), " :");
                    vector<string> poleInfo;
                    tokPoleInfo.ToVector(poleInfo);

                    int mpoleIdx = boost::lexical_cast<int>(poleInfo[0]);
                    string mpoleName = poleInfo[1];

                    mpoleIdcs.push_back(mpoleIdx);
                    mpoleNames.push_back(mpoleName);

                }

                alloc_frag_mpoleIdx[mapKeyName]         = mpoleIdcs;
                alloc_frag_mpoleName[mapKeyName]        = mpoleNames;
                alloc_frag_trihedron[mapKeyName]        = trihedron_mps;
                alloc_frag_weights[mapKeyName]          = weights_mps;
            }
        }
    }

    _alloc_frag_mpoleIdx    = alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName   = alloc_frag_mpoleName;
    _alloc_frag_trihedron   = alloc_frag_trihedron;
    _alloc_frag_weights     = alloc_frag_weights;
}


void XQMP::Collect_MPS(Topology *top) {

    // +++++++++++++ //
    // Parse + Store //
    // +++++++++++++ //

    map<int,string> ::iterator misit;

    // Neutral species
    for (misit = _segId_mpsFile_n.begin();
         misit != _segId_mpsFile_n.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, 0); }
    }

    // Anion
    for (misit = _segId_mpsFile_e.begin();
         misit != _segId_mpsFile_e.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, -1); }
    }

    // Cation
    for (misit = _segId_mpsFile_h.begin();
         misit != _segId_mpsFile_h.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    
    // Job Seg1 Seg2
    vector<XJob*> :: iterator jit;
    for (jit = _XJobs.begin();
         jit < _XJobs.end();
         ++jit) {

        string mpsfile1 = (*jit)->getMPS1();
        string mpsfile2 = (*jit)->getMPS2();

        if (_mpsFile_pSites_job.count(mpsfile1) > 0 ) { ; }
        else { _mpsFile_pSites_job[mpsfile1] = Parse_GDMA(mpsfile1, 0); }

        if (_mpsFile_pSites_job.count(mpsfile2) > 0 ) { ; }
        else { _mpsFile_pSites_job[mpsfile2] = Parse_GDMA(mpsfile2, 0); }

    }

//    map<string, vector<APolarSite*> > ::iterator it;
//    for (it = _mpsFile_pSites_job.begin();
//         it != _mpsFile_pSites_job.end();
//         ++it) {
//         cout << endl << "KEY_" << it->first << flush;
//    }
    

    cout << endl
         << "... ... ... Parsed " << _mpsFile_pSites.size()
         << " .mps files from " << _emp_file
         << ", " << _mpsFile_pSites_job.size()
         << " .mps files from " << _job_file
         << flush;
}


void XQMP::Collect_EMP(string emp_file, Topology *top) {

    std::string line;
    std::ifstream intt;
    intt.open(emp_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

            // Line sample:
            // 1022 C60 C60_n.mps C60_e.mps C60_h.mps

            int segId       = boost::lexical_cast<int>(split[0]);
            string segName  = boost::lexical_cast<string>(split[1]);

            // Some compliance checks
            assert( split.size() == 5 );                                        // <- File format correct?
            assert( top->getSegment(segId)->getName() == segName );             // <- Input matches topology?

            _segId_mpsFile_n[segId] = split[2];
            _segId_mpsFile_e[segId] = split[3];
            _segId_mpsFile_h[segId] = split[4];

        } // Exit loop over lines
    }
    else { cout << endl << "ERROR: No such file " << emp_file << endl; }

    assert( _segId_mpsFile_n.size() == top->Segments().size() );                // <- Input for all segments?

}


void XQMP::Collect_JOB(string job_file, Topology *top) {

    QMNBList &nblist = top->NBList();

    std::string line;
    std::ifstream intt;
    intt.open(job_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

// Sample line
// # JOB_ID TAG  PAIR_ID SEG1_ID SEG1_NAME SEG1_MPS SEG2_ID SEG2_NAME SEG2_MPS  (TYPE  SITE)
//   1      E_CT 3819    182     C60       c60.mps  392     DCV       dcv.mps   (site  392)

            int jobId       = boost::lexical_cast<int>(split[0]);
            string tag      = split[1];
            int pairId      = boost::lexical_cast<int>(split[2]);

            int seg1Id      = boost::lexical_cast<int>(split[3]);
            string seg1Name = split[4];
            string seg1mps  = split[5];

            int seg2Id      = boost::lexical_cast<int>(split[6]);
            string seg2Name = split[7];
            string seg2mps  = split[8];

            string job_type = "pair";
            int    energy_site_id = -1;
            if (split.size() == 11) {
                job_type = split[9];
                energy_site_id = boost::lexical_cast<int>(split[10]);
            }

            Segment *seg1   = top->getSegment(seg1Id);
            Segment *seg2   = top->getSegment(seg2Id);

            _XJobs.push_back(new XJob(jobId,  tag,    job_type, energy_site_id,
                                      pairId, seg1Id, seg2Id,   seg1mps,
                                      seg2mps, top));

            _XJobs.back()->setType(job_type, energy_site_id);

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << job_file << endl;
           throw runtime_error("Please supply input file.");           }

    cout << endl 
         << "... ... ... Registered " << _XJobs.size() << " jobs. "
         << flush;

}


void XQMP::Create_MPOLS(Topology *top) {
    
    // Warning: Direct mapping of k > 0 multipoles to MD coordinates
    bool print_huge_map2md_warning = false;
    
    // Log warning: Symmetry = 1 and k > 0 multipoles.
    map<string,bool> warned_symm_idkey;

    // +++++++++++++++++++++++++++++++++++++ //
    // Equip TOP with distributed multipoles //
    // +++++++++++++++++++++++++++++++++++++ //

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg                = *sit;
        int segId                   = seg->getId();

        bool map2md                 = _map2md[seg->getName()];

        string mps_n                = _segId_mpsFile_n[segId];
        string mps_h                = _segId_mpsFile_h[segId];
        string mps_e                = _segId_mpsFile_e[segId];

        vector<APolarSite*> pols_n   = _mpsFile_pSites[mps_n];
        vector<APolarSite*> pols_h   = _mpsFile_pSites[mps_h];
        vector<APolarSite*> pols_e   = _mpsFile_pSites[mps_e];

        // Merge polar sites
        assert(pols_n.size() == pols_h.size());
        assert(pols_n.size() == pols_e.size());

        for (int i = 0; i < pols_n.size(); i++) {
            pols_n[i]->setQs( pols_h[i]->getQs(+1), +1 );
            pols_n[i]->setPs( pols_h[i]->getPs(+1), +1 );
        }
        for (int i = 0; i < pols_n.size(); ++i) {
            pols_n[i]->setQs(pols_e[i]->getQs(-1), -1 );
            pols_n[i]->setPs(pols_e[i]->getPs(-1), -1 );
        }

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;

            // Retrieve polar-site data for this fragment
            string idkey                = frag->getName() + seg->getName()
                                        + seg->getMolecule()->getName();
            vector<int> polesInFrag     = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag  = _alloc_frag_mpoleName.at(idkey);
            vector<double> weightsInFrag= _alloc_frag_weights.at(idkey);

            if (map2md && polesInFrag.size() != frag->Atoms().size()) {
                cout << endl
                     << "ERROR: Segment " << seg->getName()
                     << " Fragment " << frag->getName()
                     << ": MAP2MD = TRUE requires same number of polar "
                     << "sites as there are atoms to perform mapping. "
                     << endl;
                throw runtime_error("Check mapping or switch to map2md = 0");
            }

            matrix rotateMP2MD;
            vec translateMP2MD;

            // Determine transformation matrices to execute mapping
            if (!map2md) {

                vector<APolarSite*> trihedron_pol;
                vector<Atom*>      trihedron_atm;

                vector<int> trihedron_ints  = _alloc_frag_trihedron.at(idkey);
                vector<int> ::iterator iit;
                for (iit = trihedron_ints.begin();
                     iit < trihedron_ints.end();
                     ++iit) {
                    trihedron_pol.push_back(pols_n[(*iit)-1]);
                }
                
                trihedron_ints  = frag->getTrihedron();
                for (iit = trihedron_ints.begin();
                     iit < trihedron_ints.end();
                     ++iit) {
                    vector< Atom* > ::iterator ait;
                    for (ait = frag->Atoms().begin();
                         ait < frag->Atoms().end();
                         ++ait) {
                        if ((*ait)->getQMId() == (*iit)) {
                            trihedron_atm.push_back(*ait);
                        }
                    }
                }


                int symmetry = trihedron_pol.size();
                assert (trihedron_pol.size() <= trihedron_atm.size() );               

                vec xMD, yMD, zMD;
                vec xQM, yQM, zQM;

                if (symmetry == 3) {
                    vec r1MD = trihedron_atm[0]->getPos();
                    vec r2MD = trihedron_atm[1]->getPos();
                    vec r3MD = trihedron_atm[2]->getPos();
                    vec r1QM = trihedron_pol[0]->getPos();
                    vec r2QM = trihedron_pol[1]->getPos();
                    vec r3QM = trihedron_pol[2]->getPos();

                    xMD = r2MD - r1MD;
                    yMD = r3MD - r1MD;
                    xQM = r2QM - r1QM;
                    yQM = r3QM - r1QM;

                    zMD = xMD ^ yMD;
                    zQM = xQM ^ yQM;

                    yMD = zMD ^ xMD;
                    yQM = zQM ^ xQM;

                    xMD = xMD.normalize();
                    yMD = yMD.normalize();
                    zMD = zMD.normalize();
                    xQM = xQM.normalize();
                    yQM = yQM.normalize();
                    zQM = zQM.normalize();
                }

                else if (symmetry == 2) {

                    vec r1MD = trihedron_atm[0]->getPos();
                    vec r2MD = trihedron_atm[1]->getPos();
                    vec r1QM = trihedron_pol[0]->getPos();
                    vec r2QM = trihedron_pol[1]->getPos();

                    xMD = r2MD - r1MD;
                    xQM = r2QM - r1QM;

                    // Normalising not necessary, but recommendable, when doing...
                    xMD = xMD.normalize();
                    xQM = xQM.normalize();

                    vec yMDtmp = vec(0,0,0);
                    vec yQMtmp = vec(0,0,0);

    // ... this: Check whether one of the components is equal or close to
    // zero. If so, this easily gives a second leg for the trihedron.
    if      ( xMD.getX()*xMD.getX() < 1e-6 ) { yMDtmp = vec(1,0,0); }
    else if ( xMD.getY()*xMD.getY() < 1e-6 ) { yMDtmp = vec(0,1,0); }
    else if ( xMD.getZ()*xMD.getZ() < 1e-6 ) { yMDtmp = vec(0,0,1); }
    if      ( xQM.getX()*xQM.getX() < 1e-6 ) { yQMtmp = vec(1,0,0); }
    else if ( xQM.getY()*xQM.getY() < 1e-6 ) { yQMtmp = vec(0,1,0); }
    else if ( xQM.getZ()*xQM.getZ() < 1e-6 ) { yQMtmp = vec(0,0,1); }

    if ( abs(yMDtmp) < 0.5 ) {
       // All components of xMD are unequal to zero => division is safe.
       // Choose vector from plane with xMD * inPlaneVec = 0:
       double tmp_x = 1.;
       double tmp_y = 1.;
       double tmp_z = 1/xMD.getZ() * (-xMD.getX()*tmp_x - xMD.getY()*tmp_y);
       yMDtmp = vec(tmp_x, tmp_y, tmp_z);
       yMDtmp.normalize();
    }
    if ( abs(yQMtmp) < 0.5 ) {
       double tmp_x = 1.;
       double tmp_y = 1.;
       double tmp_z = 1/xQM.getZ() * (-xQM.getX()*tmp_x - xQM.getY()*tmp_y);
       yQMtmp = vec(tmp_x, tmp_y, tmp_z);
       yQMtmp.normalize();
    }

                    // Now proceed as for symmetry 3
                    zMD = xMD ^ yMDtmp;
                    yMD = zMD ^ xMD;
                    zQM = xQM ^ yQMtmp;
                    yQM = zQM ^ xQM;

                    xMD.normalize();
                    yMD.normalize();
                    zMD.normalize();
                    xQM.normalize();
                    yQM.normalize();
                    zQM.normalize();
                }

                else if (symmetry == 1) {

                    if (!warned_symm_idkey[idkey]) {
                        cout << endl << "... ... ... "
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << flush;
                        warned_symm_idkey[idkey] = true;
                    }

                    xMD = vec(1,0,0);
                    yMD = vec(0,1,0);
                    zMD = vec(0,0,1);
                    xQM = vec(1,0,0);
                    yQM = vec(0,1,0);
                    zQM = vec(0,0,1);
                }

                else {
                    cout << endl
                         << "NOTE: Invalid definition of local frame in fragment "
                         << frag->getName();
                    cout << ". Assuming point particle for mapping. "
                         << endl;
                    cout << endl
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << endl;

                    xMD = vec(1,0,0);
                    yMD = vec(0,1,0);
                    zMD = vec(0,0,1);
                    xQM = vec(1,0,0);
                    yQM = vec(0,1,0);
                    zQM = vec(0,0,1);
                }

                matrix rotMD = matrix(xMD, yMD, zMD);
                matrix rotMP = matrix(xQM, yQM, zQM);
                
                rotateMP2MD = rotMD * rotMP.Transpose();


                // ++++++++++++++++++ //
                // Transform fragment //
                // ++++++++++++++++++ //

                
                vec CoMP = vec(0.,0.,0.);
                double W = 0.0;
                for (int i = 0; i < polesInFrag.size(); ++i) {

                    double weight = weightsInFrag[i];
                    
                    vec pos = pols_n[polesInFrag[i]-1]->getPos();
                    
                    CoMP += weight*pos;
                    W += weight;

                }
                CoMP /= W;

                translateMP2MD = frag->getCoMD() - CoMP;

            }            

            // Create polar sites 
            for (int i = 0; i < polesInFrag.size(); i++) {

                string name             = namesInFrag[i];
                int poleId              = polesInFrag[i];

                APolarSite *templ        = pols_n[poleId-1];
                APolarSite *newSite      = top->AddAPolarSite(name);
                newSite->ImportFrom(templ);
                seg->AddAPolarSite(newSite);
                frag->AddAPolarSite(newSite);

                // Shift + rotate
                if (!map2md) {
                    newSite->Translate(translateMP2MD);
                    newSite->Rotate(rotateMP2MD, frag->getCoMD());
                }
                else {
                    vec mdpos = frag->Atoms()[i]->getPos();
                    newSite->setPos(mdpos);
                    if (newSite->getRank() > 0) {
                        print_huge_map2md_warning = true;
                    }
                }
            }
        }
    }

    if (print_huge_map2md_warning) {
        cout << endl << endl
             << "**************************************************************"
             << "WARNING: MAP2MD = TRUE while using higher-rank multipoles can "
             << "mess up the orientation of those multipoles if the coordinate "
             << "frame used in the .mps file does not agree with the global MD "
             << "frame. If you know what you are doing - proceed ... "
             << "**************************************************************"
             << endl;
    }

    top->setIsEStatified(true);

}


vector<APolarSite*> XQMP::Map_MPols_To_Seg(vector<APolarSite*> &pols_n, Segment *seg) {

    bool print_huge_map2md_warning = false;

    vector<APolarSite*> return_pols;
    return_pols.reserve(pols_n.size());

    int segId                   = seg->getId();
    bool map2md                 = _map2md[seg->getName()];

    vector<Fragment*> ::iterator fit;
    for (fit = seg->Fragments().begin();
         fit < seg->Fragments().end();
         ++fit) {

        Fragment *frag = *fit;

        // Retrieve polar-site data for this fragment
        string idkey                = frag->getName() + seg->getName()
                                    + seg->getMolecule()->getName();
        vector<int> polesInFrag     = _alloc_frag_mpoleIdx.at(idkey);
        vector<string> namesInFrag  = _alloc_frag_mpoleName.at(idkey);
        vector<double> weightsInFrag= _alloc_frag_weights.at(idkey);

        if (map2md && polesInFrag.size() != frag->Atoms().size()) {
            cout << endl
                 << "ERROR: Segment " << seg->getName()
                 << " Fragment " << frag->getName()
                 << ": MAP2MD = TRUE requires same number of polar "
                 << "sites as there are atoms to perform mapping. "
                 << endl;
            throw runtime_error("Check mapping or switch to map2md = 0");
        }

        matrix rotateMP2MD;
        vec translateMP2MD;

        // Determine transformation matrices to execute mapping
        if (!map2md) {

            vector<APolarSite*> trihedron_pol;
            vector<Atom*>      trihedron_atm;

            vector<int> trihedron_ints  = _alloc_frag_trihedron.at(idkey);
            vector<int> ::iterator iit;
            for (iit = trihedron_ints.begin();
                 iit < trihedron_ints.end();
                 ++iit) {
                trihedron_pol.push_back(pols_n[(*iit)-1]);
            }

            trihedron_ints  = frag->getTrihedron();
            for (iit = trihedron_ints.begin();
                 iit < trihedron_ints.end();
                 ++iit) {
                vector< Atom* > ::iterator ait;
                for (ait = frag->Atoms().begin();
                     ait < frag->Atoms().end();
                     ++ait) {
                    if ((*ait)->getQMId() == (*iit)) {
                        trihedron_atm.push_back(*ait);
                    }
                }
            }


            int symmetry = trihedron_pol.size();
            assert (trihedron_pol.size() == trihedron_atm.size() );

            vec xMD, yMD, zMD;
            vec xQM, yQM, zQM;

            if (symmetry == 3) {
                vec r1MD = trihedron_atm[0]->getPos();
                vec r2MD = trihedron_atm[1]->getPos();
                vec r3MD = trihedron_atm[2]->getPos();
                vec r1QM = trihedron_pol[0]->getPos();
                vec r2QM = trihedron_pol[1]->getPos();
                vec r3QM = trihedron_pol[2]->getPos();

                xMD = r2MD - r1MD;
                yMD = r3MD - r1MD;
                xQM = r2QM - r1QM;
                yQM = r3QM - r1QM;

                zMD = xMD ^ yMD;
                zQM = xQM ^ yQM;

                yMD = zMD ^ xMD;
                yQM = zQM ^ xQM;

                xMD = xMD.normalize();
                yMD = yMD.normalize();
                zMD = zMD.normalize();
                xQM = xQM.normalize();
                yQM = yQM.normalize();
                zQM = zQM.normalize();
            }

            else if (symmetry == 2) {

                vec r1MD = trihedron_atm[0]->getPos();
                vec r2MD = trihedron_atm[1]->getPos();
                vec r1QM = trihedron_pol[0]->getPos();
                vec r2QM = trihedron_pol[1]->getPos();

                xMD = r2MD - r1MD;
                xQM = r2QM - r1QM;

                // Normalising not necessary, but recommendable, when doing...
                xMD = xMD.normalize();
                xQM = xQM.normalize();

                vec yMDtmp = vec(0,0,0);
                vec yQMtmp = vec(0,0,0);

    // ... this: Check whether one of the components is equal or close to
    // zero. If so, this easily gives a second leg for the trihedron.
    if      ( xMD.getX()*xMD.getX() < 1e-6 ) { yMDtmp = vec(1,0,0); }
    else if ( xMD.getY()*xMD.getY() < 1e-6 ) { yMDtmp = vec(0,1,0); }
    else if ( xMD.getZ()*xMD.getZ() < 1e-6 ) { yMDtmp = vec(0,0,1); }
    if      ( xQM.getX()*xQM.getX() < 1e-6 ) { yQMtmp = vec(1,0,0); }
    else if ( xQM.getY()*xQM.getY() < 1e-6 ) { yQMtmp = vec(0,1,0); }
    else if ( xQM.getZ()*xQM.getZ() < 1e-6 ) { yQMtmp = vec(0,0,1); }

    if ( abs(yMDtmp) < 0.5 ) {
    // All components of xMD are unequal to zero => division is safe.
    // Choose vector from plane with xMD * inPlaneVec = 0:
    double tmp_x = 1.;
    double tmp_y = 1.;
    double tmp_z = 1/xMD.getZ() * (-xMD.getX()*tmp_x - xMD.getY()*tmp_y);
    yMDtmp = vec(tmp_x, tmp_y, tmp_z);
    yMDtmp.normalize();
    }
    if ( abs(yQMtmp) < 0.5 ) {
    double tmp_x = 1.;
    double tmp_y = 1.;
    double tmp_z = 1/xQM.getZ() * (-xQM.getX()*tmp_x - xQM.getY()*tmp_y);
    yQMtmp = vec(tmp_x, tmp_y, tmp_z);
    yQMtmp.normalize();
    }

                // Now proceed as for symmetry 3
                zMD = xMD ^ yMDtmp;
                yMD = zMD ^ xMD;
                zQM = xQM ^ yQMtmp;
                yQM = zQM ^ xQM;

                xMD.normalize();
                yMD.normalize();
                zMD.normalize();
                xQM.normalize();
                yQM.normalize();
                zQM.normalize();
            }

            else if (symmetry == 1) {

                //cout << endl
                //     << "WARNING: Symmetry = 1 for fragment "
                //     << frag->getName() << ": This will generate artifacts "
                //     << "when mapping higher-rank multipoles (dipoles, ..)."
                //     << endl;

                xMD = vec(1,0,0);
                yMD = vec(0,1,0);
                zMD = vec(0,0,1);
                xQM = vec(1,0,0);
                yQM = vec(0,1,0);
                zQM = vec(0,0,1);
            }

            else {
                cout << endl
                     << "NOTE: Invalid definition of local frame in fragment "
                     << frag->getName();
                cout << ". Assuming point particle for mapping. "
                     << endl;
                cout << endl
                     << "WARNING: Symmetry = 1 for fragment "
                     << frag->getName() << ": This will generate artifacts "
                     << "when mapping higher-rank multipoles (dipoles, ..)."
                     << endl;

                xMD = vec(1,0,0);
                yMD = vec(0,1,0);
                zMD = vec(0,0,1);
                xQM = vec(1,0,0);
                yQM = vec(0,1,0);
                zQM = vec(0,0,1);
            }

            matrix rotMD = matrix(xMD, yMD, zMD);
            matrix rotMP = matrix(xQM, yQM, zQM);

            rotateMP2MD = rotMD * rotMP.Transpose();


            // ++++++++++++++++++ //
            // Transform fragment //
            // ++++++++++++++++++ //


            vec CoMP = vec(0.,0.,0.);
            double W = 0.0;
            for (int i = 0; i < polesInFrag.size(); ++i) {

                double weight = weightsInFrag[i];

                vec pos = pols_n[polesInFrag[i]-1]->getPos();

                CoMP += weight*pos;
                W += weight;

            }
            CoMP /= W;

            translateMP2MD = frag->getCoMD() - CoMP;

        }

        // Create polar sites
        for (int i = 0; i < polesInFrag.size(); i++) {

            string name             = namesInFrag[i];
            int poleId              = polesInFrag[i];

            APolarSite *templ        = pols_n[poleId-1];
            APolarSite *newSite      = new APolarSite(-1, name);
            newSite->ImportFrom(templ);

            // Shift + rotate
            if (!map2md) {
                newSite->Translate(translateMP2MD);
                newSite->Rotate(rotateMP2MD, frag->getCoMD());
            }
            else {
                vec mdpos = frag->Atoms()[i]->getPos();
                newSite->setPos(mdpos);
                if (newSite->getRank() > 0) {
                    print_huge_map2md_warning = true;
                }
            }

            newSite->Charge(0);
            
            return_pols.push_back(newSite);

        }
    } // End loop over fragments

    if (print_huge_map2md_warning) {
        cout << endl << endl
             << "**************************************************************"
             << "WARNING: MAP2MD = TRUE while using higher-rank multipoles can "
             << "mess up the orientation of those multipoles if the coordinate "
             << "frame used in the .mps file does not agree with the global MD "
             << "frame. If you know what you are doing - proceed ... "
             << "**************************************************************"
             << endl;
    }

    return return_pols;
}


vector<APolarSite*> XQMP::Parse_GDMA(string filename, int state) {

    return APS_FROM_MPS(filename,state);    
    
}


bool XQMP::EvaluateFrame(Topology *top) {

    cout << endl
         << "... ... Load multipole definition, collect jobs: "
         << flush;

    Collect_XML(_xml_file, top); // <- Allocate polar sites to fragments    
    Collect_EMP(_emp_file, top); // <- Collect segs .mps => Background .mps
    Collect_JOB(_job_file, top); // <- Collect jobs .mps => Foreground .mps
    Collect_MPS(top);            // <- Parse all .mps files (fore-, background)


    // ++++++++++++++++++++ //
    // Ridigidfy + Estatify //
    // ++++++++++++++++++++ //

    // Rigidify segments, fragments
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else {
        cout << endl
             << "... ... System is already rigidified.";
    }

    // Create polar sites
    if (top->isEStatified() == false) {
        this->Create_MPOLS(top);
        cout << endl
             << "... ... Created " << top->APolarSites().size()
             << " multipole sites. "
             << flush;
    }
    else {
        cout << endl
             << "... ... System is already estatified. "
             << "This prohibits running XQMultipole. "
             << flush;
    }

    // To check rotation into global frame
    if (_pdb_check != "") {
        cout << endl
             << "... ... Writing polar-site coordinates to "
             << _pdb_check << ". "
             << flush;

        string mpNAME   = _pdb_check;
        FILE *mpPDB     = NULL;
        mpPDB           = fopen(mpNAME.c_str(), "w");

        vector<Segment*>::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             ++sit) {
            (*sit)->WritePDB(mpPDB, "Multipoles", "Charges");
        }
        fclose(mpPDB);
    }

    
    // +++++++++++++++++++++++++++++++++ //
    // Create + start threads (Job Ops)  //
    // +++++++++++++++++++++++++++++++++ //

    // Convert threads into subthreads if beneficial
    if (_XJobs.size() < _nThreads) {
        _subthreads = (_nThreads - _XJobs.size()) / _XJobs.size() + 1;
        _nThreads   = _XJobs.size();

        cout << endl << "... ... "
             << "Converted threads into subthreads to increase efficiency: "
             << "NT = " << _nThreads << ", NST = " << _subthreads
             << flush;
    }

    vector<JobXQMP*> jobOps;    

    _nextXJob = _XJobs.begin();

    for (int id = 0; id < _nThreads; id++) {
        JobXQMP *newOp = new JobXQMP(id, top, this);
        jobOps.push_back(newOp);
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->InitSlotData(top);
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->Start();
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->WaitDone();
    }

    for (int id = 0; id < _nThreads; id++) {
        delete jobOps[id];
    }

    jobOps.clear();


    FILE *out;
    out = fopen(this->_outFile.c_str(), "w");
    vector<XJob*> :: iterator jit;
    for (jit = _XJobs.begin(); jit < _XJobs.end(); ++jit) {
        (*jit)->WriteInfoLine(out);
    }
    fclose(out);

}



// ========================================================================== //
//                            JOBXQMP MEMBER FUNCTIONS                         //
// ========================================================================== //


void XQMP::JobXQMP::InitSlotData(Topology *top) {

    vector< Segment* > ::iterator sitRef;
    vector< vector<APolarSite*> > ::iterator sitNew;
    vector< APolarSite* > ::iterator pitRef;
    vector< APolarSite* > ::iterator pitNew;

    _polarSites.resize(top->Segments().size());
    assert(top->Segments().size() == _polarSites.size());

    for (sitRef = top->Segments().begin(), sitNew = _polarSites.begin();
         sitRef < top->Segments().end();
         ++sitRef, ++sitNew) {

        (*sitNew).resize((*sitRef)->APolarSites().size());

        for (pitRef = (*sitRef)->APolarSites().begin(),
             pitNew = (*sitNew).begin();
             pitRef < (*sitRef)->APolarSites().end();
             ++pitRef, ++ pitNew) {

            *pitNew = new APolarSite();
            (*pitNew)->ImportFrom(*pitRef, "full");
            (*pitNew)->Charge(0);
        }
    }
}


void XQMP::JobXQMP::Run(void) {

    while (true) {
        _job = _master->RequestNextJob(_id, _top);

        if (_job == NULL) { break; }
        else { this->EvalJob(_top, _job); }
    }
}


void XQMP::JobXQMP::EvalJob(Topology *top, XJob *job) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    // ++++++++++++++++++++++++++ //
    // Adapt polar sites          //
    // ++++++++++++++++++++++++++ //
    _polarSites_job = _polarSites;

    int subs_here1 = job->getSeg1Id();
    int subs_here2 = job->getSeg2Id();

    vector<APolarSite*> subs1_raw = _master->_mpsFile_pSites_job[job->getMPS1()];
    vector<APolarSite*> subs2_raw = _master->_mpsFile_pSites_job[job->getMPS2()];
    vector<APolarSite*> subs1 = _master->Map_MPols_To_Seg(subs1_raw,job->Seg1());
    vector<APolarSite*> subs2 = _master->Map_MPols_To_Seg(subs2_raw,job->Seg2());

    _polarSites_job[subs_here1-1] = subs1;
    _polarSites_job[subs_here2-1] = subs2;
    
    

    // ++++++++++++++++++++++++++ //
    // Define polarization sphere //
    // ++++++++++++++++++++++++++ //

    vec center = job->Center();
    // vec center = job->Seg1()->getPos(); // UNCOMMENT TO COMPARE TO EMULTIPOLE

    this->_segsPolSphere.clear(); // <- Segments    within cutoff
    this->_segsOutSphere.clear(); // <- Segments    within cutoff1, cutoff2
    this->_polsPolSphere.clear(); // <- Polar sites within cutoff
    this->_polsOutSphere.clear(); // <- Polar sites within cutoff1, cutoff2

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

        double r12 = abs(_top->PbShortestConnect((*sit)->getPos(), center));

        // Always add pair-job segments to polSphere, even for cut-off = 0.0
        if ( (*sit)->getId() == job->getSeg1Id()
          || (*sit)->getId() == job->getSeg2Id()) {
            if   (job->getType() == "pair") { r12 = -1; }
            else                            { ; }
        }

        if      ( r12 > _master->_cutoff2) { continue; }

        else if ( r12 > _master->_cutoff ) {
            _segsOutSphere.push_back(*sit);
            _polsOutSphere.push_back( _polarSites_job[(*sit)->getId() - 1] );
        }
        else {
            _segsPolSphere.push_back(*sit);
            _polsPolSphere.push_back( _polarSites_job[(*sit)->getId() - 1] );
        }
    }

    if (_master->_maverick) {
        cout << endl
             << "... ... ... Segments in polarization sphere: "
             << _segsPolSphere.size()
             << "; segments in static shell: "
             << _segsOutSphere.size()
             << flush;
    }

    job->setSizePol(_polsPolSphere.size());
    job->setSizeShell(_polsOutSphere.size());

//    FILE *out;
//    string shellFile = "OuterShell.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    for (sit = _segsOutSphere.begin(); sit < _segsOutSphere.end(); ++sit) {
//        (*sit)->WritePDB(out, "Multipoles", "");
//    }
//    fclose(out);
//
//    shellFile = "InnerShell.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    for (sit = _segsPolSphere.begin(); sit < _segsPolSphere.end(); ++sit) {
//        (*sit)->WritePDB(out, "Multipoles", "");
//    }
//    fclose(out);
//
//    shellFile = "Pair.pdb";
//    out = fopen(shellFile.c_str(), "w");
//    job->Seg1()->WritePDB(out, "Multipoles", "");
//    job->Seg2()->WritePDB(out, "Multipoles", "");
//    fclose(out);


    _inductor = XInductor(_master->_induce,     _master->_induce_intra_pair,
                          _master->_subthreads, _master->_wSOR_N,
                          _master->_epsTol,     _master->_maxIter,
                          _master->_maverick,   top,     _id);
    
    _inductor.Evaluate(job, _segsPolSphere, _segsOutSphere,
                            _polsPolSphere, _polsOutSphere,
                            _polarSites,    _polarSites_job);
    

    // ++++++++++++++++++++++++++ //
    // Write Checkpoint File      //
    // ++++++++++++++++++++++++++ //

    if (_master->_write_chk) {

        string dotsuffix = "";

        if (_master->_chk_format == "gaussian") {
            dotsuffix = ".com";
        }
        else if (_master->_chk_format == "xyz") {
            dotsuffix = ".xyz";
        }

        FILE *out;
        string chk_file = job->getTag()+_master->_write_chk_suffix+dotsuffix;
        out = fopen(chk_file.c_str(), "w");

        vector< APolarSite* >         ::iterator pit;
        int pcount = 0;

        // Count polar sites for header line in xyz
        if (_master->_chk_format == "xyz") {
            for (sit = _segsPolSphere.begin(); 
                 sit < _segsPolSphere.end();
                 ++sit) {

                int segId = (*sit)->getId();

                if (job->isInCenter(segId)) {
                    pcount += _polarSites_job[segId-1].size();
                }
                else {
                    pcount += 3*_polarSites_job[segId-1].size();
                }
            }
            for (sit = _segsOutSphere.begin(); 
                 sit < _segsOutSphere.end();
                 ++sit) {

                int segId = (*sit)->getId();

                if (job->isInCenter(segId)) {
                    assert(false); // Central region in outer shell!? Error.
                }
                else {
                    pcount += _polarSites_job[segId-1].size();
                }
            }
            
            fprintf(out, "%1d\n"
                         "XYZ WITH DIPOLES IN POLARIZABLE "
                         "SHELL SPLIT ONTO POINT CHARGES\n", pcount);
        }

        
        // Save coordinates of central region            
        for (sit = _segsPolSphere.begin(); 
             sit < _segsPolSphere.end(); ++sit) {

            int segId = (*sit)->getId();

            if (job->isInCenter(segId)) {

                vec pb_shift = job->Center() - (*sit)->getPos()
                     - top->PbShortestConnect((*sit)->getPos(), job->Center());

                for (pit = _polarSites_job[segId-1].begin();
                     pit < _polarSites_job[segId-1].end();
                     ++pit) {
                    (*pit)->WriteXyzLine(out, pb_shift, _master->_chk_format);
                }
            }
        }

        if (_master->_chk_format == "gaussian") {
            fprintf(out, "\n");
        }

        // Save induction state of polarizable sphere
        for (sit = _segsPolSphere.begin(); 
             sit < _segsPolSphere.end(); ++sit) {

            int segId = (*sit)->getId();

            if (job->isInCenter(segId)) { continue; }

            vec pb_shift = job->Center() - (*sit)->getPos()
                     - top->PbShortestConnect((*sit)->getPos(), job->Center());

            for (pit = _polarSites_job[segId-1].begin();
                 pit < _polarSites_job[segId-1].end();
                 ++pit) {
                 (*pit)->WriteChkLine(out, pb_shift,
                                           _master->_chk_split_dpl,
                                           _master->_chk_format,
                                           _master->_chk_dpl_spacing);
            }
        }

        // Write point charges of outer sphere
        for (sit = _segsOutSphere.begin();
             sit < _segsOutSphere.end(); ++sit) {

            int segId = (*sit)->getId();

            if (job->isInCenter(segId)) {
                    assert(false); // Central region in outer shell!? Error.             
            }

            vec pb_shift = job->Center() - (*sit)->getPos()
                     - top->PbShortestConnect((*sit)->getPos(), job->Center());

            for (pit = _polarSites_job[segId-1].begin();
                 pit < _polarSites_job[segId-1].end();
                 ++pit) {
                 (*pit)->WriteChkLine(out, pb_shift,
                                           false,
                                           _master->_chk_format,
                                           _master->_chk_dpl_spacing);
            }
        }

        fclose(out);
    }

    // ++++++++++++++++++++++++++ //
    // Clean up polar sites       //
    // ++++++++++++++++++++++++++ //

    vector< APolarSite* > ::iterator cleanit;
    for (cleanit = subs1.begin(); cleanit < subs1.end(); ++cleanit) {
        delete *cleanit;
    }
    for (cleanit = subs2.begin(); cleanit < subs2.end(); ++cleanit) {
        delete *cleanit;
    }
    subs1.clear();
    subs2.clear();

}


XJob *XQMP::RequestNextJob(int id, Topology *top) {

    _nextJobMutex.Lock();

    XJob *workOnThis;

    if (_nextXJob == _XJobs.end()) {
        workOnThis = NULL;
    }
    else {
        workOnThis = *_nextXJob;
        _nextXJob++;
        cout << endl 
             << "... ... Thread " << id << " evaluating job "
             << workOnThis->getId() << " " << workOnThis->getTag()
             << flush;
    }

    _nextJobMutex.Unlock();

    return workOnThis;
}







}}

#endif

