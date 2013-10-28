#ifndef VOTCA_CTP_LOG2MPS_H
#define VOTCA_CTP_LOG2MPS_H


#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>
#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/qmmachine.h>


namespace votca { namespace ctp {


class Log2Mps : public QMTool
{
public:

    Log2Mps() { };
   ~Log2Mps() { };

    string Identify() { return "log2mps"; }

    void   Initialize(Property *options);
    bool   Evaluate();



private:

    string _package;
    string _logfile;
    string _mpsfile;

};


void Log2Mps::Initialize(Property *opt) {
    
    QMPackageFactory::RegisterAll();
    
    string key = "options.log2mps";
    _package = opt->get(key+".package").as<string>();
    _logfile = opt->get(key+".logfile").as<string>();
    

    _mpsfile = (opt->exists(key+".mpsfile")) ? 
        opt->get(key+".mpsfile").as<string>() : "";
    if (_mpsfile == "") _mpsfile = _logfile.substr(0,_logfile.size()-4)+".mps";

    cout << endl << "... ... " << _logfile << " => " << _mpsfile << flush;
}


bool Log2Mps::Evaluate() {
    
    // Logger (required for QM package, so we can just as well use it)
    Logger log;
    log.setPreface(logINFO, "\n... ...");
    log.setPreface(logDEBUG, "\n... ...");
    log.setReportLevel(logDEBUG);
    log.setMultithreading(true);  
    
    // Set-up QM package
    LOG(logINFO,log) << "Using package <" << _package << ">" << flush;
    QMPackage *qmpack = QMPackages().Create(_package);    
    qmpack->doGetCharges(true);
    qmpack->setLog(&log);
    qmpack->setRunDir(".");
    qmpack->setLogFileName(_logfile);
    
    // Create orbitals, fill with life & extract QM atoms
    Orbitals orbs;
    int cdx = qmpack->ParseLogFile(&orbs);
    if (!cdx) {
        cout << "\nERROR Parsing " << _logfile << "failed. Abort." << endl;
        throw std::runtime_error("(see above, parsing error)");
    }    
    vector<QMAtom*> &qmatoms = *orbs.getAtoms();
    vector<QMAtom*>::iterator it;
    
    // Sanity checks, total charge
    double Q = 0.0;
    for (it = qmatoms.begin(); it < qmatoms.end(); ++it) {
        Q += (*it)->charge;
    }
    
    if (qmatoms.size() < 1) {
        cout << "\nERROR No charges extracted from " << _logfile 
            << ". Abort.\n" << flush;
        throw std::runtime_error("(see above, input or parsing error)");
    }
    LOG(logINFO,log) 
        << qmatoms.size() << " QM atoms, total charge Q = " << Q << flush;    
    
    
    // Convert to polar segment & write mps-file
    QMMInterface qmmface;
    PolarSeg *pseg = qmmface.Convert(qmatoms);
    
    string tag = "::LOG2MPS " 
        + (boost::format("(log-file='%1$s' : %2$d QM atoms)")
        % _logfile % qmatoms.size()).str();    
    pseg->WriteMPS(_mpsfile, tag);
}



}}

#endif