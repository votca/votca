// Local private VOTCA includes
#include "ewald.h"


namespace votca {
namespace xtp {


template<class EwaldMethod>
void Ewald<EwaldMethod>::Initialize(tools::Property *opt) {

    // NOTE These options are passed on to <EwaldMethod> in ::EvalJob
    _options = opt;
    //_maverick = (_nThreads == 1) ? true : false;
    
    std::cout << std::endl
         << "... ... Initialized with " << nThreads_ << " threads. "
         << std::endl;

    std::string key = "options.ewald.jobcontrol";
        if ( opt->exists(key+".job_file")) {
            jobfile_ = opt->get(key+".job_file").as<std::string>();
        }
        else {
            std::cout << std::endl;
            throw std::runtime_error("Job-file not set. Abort.");
        }    
    
    key = "options.ewald.multipoles";
        if (opt->exists(key+".mapping")) {
            _xml_file = opt->get(key+".mapping").as< std::string >();
        }
        else {
            std::cout << std::endl;
            throw std::runtime_error("Multipole mapping file not set. Abort.");
        }
        if ( opt->exists(key+".mps_table")) {
            _mps_table = opt->get(key+".mps_table").as<std::string>();
        }
        else {
            std::cout << std::endl;
            throw std::runtime_error("Background mps table not set. Abort.");
        }
        if ( opt->exists(key+".polar_bg")) {
            _polar_bg_arch = opt->get(key+".polar_bg").as<std::string>();
        }
        else { _polar_bg_arch = ""; }
        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<bool>();
        }
        else { _pdb_check = false; }
        if (opt->exists(key+".ptop_check")) {
            _ptop_check = opt->get(key+".ptop_check").as<bool>();
        }
        else { _ptop_check = false; }
    
    return;
}


template<class EwaldMethod>
void Ewald<EwaldMethod>::PreProcess(Topology *top) {
    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    std::cout << std::endl << "... ... Initialize MPS-mapper: " << std::flush;
    _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
    return;
}


template<class EwaldMethod>
void Ewald<EwaldMethod>::ReadJobFile(Topology& top) {    
    assert(false && "<::ReadJobFile> NOT IMPLEMENTED");    
    return;
}


template<class EwaldMethod>
void Ewald<EwaldMethod>::WriteJobFile(const Topology& top) {
    
    
    // SET UP FILE STREAM
    std::ofstream ofs;
    std::string jobFile = "ewald_jobs.xml";
    ofs.open(jobFile.c_str(), std::ofstream::out);
    if (!ofs.is_open()) throw std::runtime_error("Bad file handle: " + jobFile);
    
    ofs << "<jobs>" << std::endl;
    
    int jobCount = 0;    
    //std::vector<Segment*>::iterator sit1;
    
    // DEFINE PAIR CHARGE STATES
    std::vector<std::string > states;
    std::vector<std::string> ::iterator vit;
    states.push_back("n");
    states.push_back("e");
    states.push_back("h");
    
    // CREATE JOBS FOR ALL SEGMENTS AND STATES
    //for (sit1 = top.Segments().begin(); sit1 < top.Segments().end(); ++sit1) {
    for (auto& seg1 : top.Segments()){// sit1 = top.Segments().begin(); sit1 < top.Segments().end(); ++sit1) {
        //Segment *seg1 = *sit1;

        int id1 = seg1.getId();
        std::string name1 = seg1.getType(); 
        
        for (vit = states.begin(); vit != states.end(); ++vit) {
            int id = ++jobCount;
            std::string s1 = *vit;
            std::string tag = (format("%1$d:%3$s:%2$s") % id1 % s1 % name1).str();
            
            tools::Property input;
            tools::Property &out = input.add("input","");
            tools::Property *next = NULL;
            next = &out.add("segment", "");
            next->add("id", (format("%1$d") % id1).str());
            next->add("type", (format("%1$s") % name1).str());
            next->add("mps", (format("MP_FILES/%1$s_%2$s.mps") 
                % name1 % s1).str());
            
            Job job(id, tag, input, Job::AVAILABLE);
            job.ToStream(ofs);
        }
    }
    
    // CLOSE STREAM
    ofs << "</jobs>" << std::endl;    
    ofs.close();
}


template<class EwaldMethod>
XJob Ewald<EwaldMethod>::ProcessInputString(Job &job, Topology &top, 
    QMThread &thread) {
    
    // Input std::string looks like this:
    // <id1>:<name1>:<mpsfile1> <id2>:<name2>: ... ... ...
    
    std::string input = job.getInput().as<std::string>();
    std::vector<Segment*> qmSegs;
    std::vector<std::string>   qmSegMps;
    std::vector<std::string> split;
    tools::Tokenizer toker(input, " \t\n");
    split = toker.ToVector();

    for (unsigned int i = 0; i < split.size(); ++i) {
                
        std::string id_seg_mps = split[i];
        std::vector<std::string> split_id_seg_mps;
        tools::Tokenizer toker(id_seg_mps, ":");
        split_id_seg_mps = toker.ToVector();

        int segId = boost::lexical_cast<int>(split_id_seg_mps[0]);
        std::string segName = split_id_seg_mps[1];
        std::string mpsFile = split_id_seg_mps[2];

        Segment seg = top.getSegment(segId);
        if (seg.getType() != segName) {
            XTP_LOG(Log::error,thread.getLogger())
                << "ERROR: Seg " << segId << ":" << seg.getType() << " "
                << " maltagged as " << segName << ". Skip job ..." << std::flush;
            throw std::runtime_error("Input does not match topology.");
        }

        qmSegs.push_back(&seg);
        qmSegMps.push_back(mpsFile);               
    }
    
    return XJob(job.getId(), job.getTag(), qmSegs, qmSegMps, &top);
}


template<class EwaldMethod>
Job::JobResult Ewald<EwaldMethod>::EvalJob(const xtp/include/votca/xtp/ewald/ewaldnd.hTopology& top, Job& job,
    QMThread& thread) {
    
    boost::timer::cpu_timer cpu_t;
    cpu_t.start();
    boost::timer::cpu_times t_in = cpu_t.elapsed();
    
    Logger *log = thread->getLogger();    
    CTP_LOG(logINFO,*log)
        << "Job input = " << job->getInput().as<std::string>() << std::flush;
    
    // CREATE XJOB FROM JOB INPUT std::string
    XJob xjob = this->ProcessInputstd::string(job, top, thread);    
    
    // GENERATE POLAR TOPOLOGY (GENERATE VS LOAD IF PREPOLARIZED)
    if (_polar_bg_arch == "") {
        CTP_LOG(logINFO,*log) << "Mps-Mapper: Generate FGC FGN BGN" << std::flush;
        _mps_mapper.Gen_FGC_FGN_BGN(top, &xjob, thread);
    }
    else {
        CTP_LOG(logINFO,*log) << "Mps-Mapper: Generate FGC, load FGN BGN from '" 
                << _polar_bg_arch << "'" << std::flush;
        _mps_mapper.Gen_FGC_Load_FGN_BGN(top, &xjob, _polar_bg_arch, thread);
    }
    
    // CALL THOLEWALD MAGIC
    EwaldMethod ewaldnd = EwaldMethod(top, xjob.getPolarTop(), _options, 
        thread->getLogger());
    if (_pdb_check)
        ewaldnd.WriteDensitiesPDB(xjob.getTag()+".densities.pdb");
    ewaldnd.Evaluate();
    if (_ptop_check)
        ewaldnd.WriteDensitiesPtop(xjob.getTag()+".fg.ptop", 
            xjob.getTag()+".mg.ptop", xjob.getTag()+".bg.ptop");
    
    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    Property output = ewaldnd.GenerateOutputstd::string();
    Job::JobResult jres = Job::JobResult();
    jres.setOutput(output);
    jres.setStatus(Job::COMPLETE);
    
    if (!ewaldnd.Converged()) {
        jres.setStatus(Job::FAILED);
        jres.setError(ewaldnd.GenerateErrorstd::string());
        CTP_LOG(logERROR,*log) << ewaldnd.GenerateErrorstd::string() << std::flush;
    }
    
    boost::timer::cpu_times t_out = cpu_t.elapsed();
    double t_run = (t_out.wall-t_in.wall)/1e9/60.;
    CTP_LOG(logINFO,*log)
        << "Job runtime was " << t_run << " min" << std::flush;
    
    return jres;
}

}}