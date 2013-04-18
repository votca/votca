#ifndef XQMULTIPOLE_H
#define XQMULTIPOLE_H


#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/xinteractor.h>



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
    XMpsMap                        _mps_mapper;
    
    // Job info : JOB_ID PAIR_ID MPS_1 MPS_2 TAG
    vector<XJob*>                  _XJobs;
    vector<XJob*>::iterator        _nextXJob;
    Mutex                          _nextJobMutex;
    Mutex                          _coutMutex;
    Mutex                          _logMutex;
    bool                           _maverick;

    /*
    vector<int>                    _jobIds;
    map<int,string>                _jobId_jobTag;
    map<int,int>                   _jobId_pairId;
    map<int,string>                _jobId_mpsFile1;
    map<int,string>                _jobId_mpsFile2;
    map<int, pair<int,int> >       _pairId_seg1Id_seg2Id;
    vector<int>::iterator          _nextJob;
    */
    
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


bool XQMP::EvaluateFrame(Topology *top) {

    cout << endl
         << "... ... Load multipole definition, collect jobs: "
         << flush;

    Collect_JOB(_job_file, top); // <- Collect jobs .mps => Foreground .mps
    
    _mps_mapper.CollectMapFromXML(_xml_file);
    _mps_mapper.CollectSegMpsAlloc(_emp_file, top);    
    _mps_mapper.CollectSitesFromMps(_XJobs);

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
        _mps_mapper.EquipWithPolSites(top);
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

    vector<APolarSite*> subs1_raw = _master->_mps_mapper.GetRawPolSitesJob(job->getMPS1());
    vector<APolarSite*> subs2_raw = _master->_mps_mapper.GetRawPolSitesJob(job->getMPS2());
    vector<APolarSite*> subs1 = _master->_mps_mapper.MapPolSitesToSeg(subs1_raw,job->Seg1());
    vector<APolarSite*> subs2 = _master->_mps_mapper.MapPolSitesToSeg(subs2_raw,job->Seg2());

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

