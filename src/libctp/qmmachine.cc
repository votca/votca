#include <votca/ctp/qmmachine.h>
#include <sys/stat.h>


namespace votca { namespace ctp {


QMMachine::QMMachine(XJob *job, XInductor *xind, Gaussian *qmpack,
                     Property *opt, string sfx, int nst, bool mav) 
                   : _job(job), _xind(xind), _qmpack(qmpack), _subthreads(nst) {
    
    string key = sfx + ".qmmm";
    
}


void QMMachine::Evaluate(XJob *job) {
    
    cout << endl << "... ... ... QM-MM Machine started." << flush;
    
    // PREPARE JOB DIRECTORY
    string jobFolder = "xjob_" + boost::lexical_cast<string>(_job->getId())
                     + "_" + _job->getTag();    
    mkdir(jobFolder.c_str(), 0755);
    
    // FIGURE OUT CHARGE + MULTIPLICITY
    double dQ = 0.0;
    for (int i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
        dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
    }
    int chrg = round(dQ);
    int spin = ( (chrg < 0) ? -chrg:chrg ) % 2 + 1;
    cout << endl << "... ... ... Q = " << chrg << ", 2S+1 = " << spin << flush;
    
    // SET ITERATION-TIME CONSTANTS
    _qmpack->setCharge(chrg);
    _qmpack->setSpin(spin);
    _qmpack->setThreads(_subthreads);
    
    // ITERATIONS
    int iter = 0;
    string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);
    mkdir(runFolder.c_str(), 0755);    
    
    // RUN CLASSICAL INDUCTION
    _xind->Evaluate(_job);
    
    // WRITE AND SET INPUT FILE, RUN
    string comFile = _job->getTag() + ".com";
    string logFile = _job->getTag() + ".log";
    string path_comFile = runFolder + "/" + comFile;
    string path_logFile = runFolder + "/" + logFile;
    this->WriteQMPackInputFile(path_comFile, _qmpack, _job);
    _qmpack->setRunDir(runFolder);
    _qmpack->setInputFile(comFile);
    
    // RUN HERE (OVERRIDE - COPY EXISTING LOG-FILE)
    string cpstr = "cp e_1_n.log " + path_logFile;
    int sig = system(cpstr.c_str());
    _qmpack->setLogFile(path_logFile);
    
    // EXTRACT LOG-FILE INFOS
    Orbitals orb_iter;
    _qmpack->ParseLogFile(&orb_iter);
    
    assert(orb_iter.hasSelfEnergy());
    assert(orb_iter.hasQMEnergy());
    
    cout << endl << "... ... ... SF energy = " << orb_iter.getSelfEnergy() << flush;
    cout << endl << "... ... ... QM energy = " << orb_iter.getQMEnergy() << flush;
    
    // while not converged:
    // ... run inductor
    // ... convert into com-file
    // ... run gaussian
    // ... read logfile
    // ... assert convergence
    // ... modify polar topology
    
    
    // mkdir(runDir) 
    // qmp setRunDir setInputFile
    // qmp WriteInputFile Run
    // qmp setLogFile ParseLogFile    
    
    return;
}


void QMMachine::WriteQMPackInputFile(string inputFile, Gaussian *qmpack, XJob *job) {
    
    // INPUT FILE
    FILE *out;
    out = fopen(inputFile.c_str(), "w");
    _qmpack->WriteInputHeader(out, job->getTag());
    job->getPolarTop()->PrintInduState(out, _qmpack->getPackageName(), true, 1e-04);
    fclose(out);
    
}


void QMMachine::ConvertPSitesToQMAtoms(vector< PolarSeg* > &psegs,
                                       vector< QMAtom * > &qmatoms) {
    
    assert(qmatoms.size() == 0);
    
    return;   
}


void QMMachine::ConvertQMAtomsToPSites(vector< QMAtom* > &qmatoms,
                                       vector< PolarSeg* > &psegs) {
    return;
}


   
    
    
    
    
    
}}
