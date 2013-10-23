#ifndef VOTCA_CTP_LOG2MPS_H
#define VOTCA_CTP_LOG2MPS_H


#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>
#include <votca/ctp/qmpackagefactory.h>

//#include <votca/ctp/orbitals.h>


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
    
}


bool Log2Mps::Evaluate() {
    
    cout << endl << "Package <" << _package << "> : " 
        << _logfile << " => " << _mpsfile << flush;
    
    
    
    Orbitals orbs;
    
    QMPackage *qmpack = QMPackages().Create(_package);
    Logger log;
    qmpack->setLog(&log);
    int cdx = qmpack->ParseLogFile(&orbs);
    
    if (!cdx) {
        cout << endl << "ERROR Parsing " << _logfile << "failed. Abort." << endl;
        throw std::runtime_error("(see above, parsing error)");
    }
    

    
}



}}

#endif