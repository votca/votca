/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef VOTCA_XTP_LOG2MPS_H
#define VOTCA_XTP_LOG2MPS_H


#include <boost/format.hpp>
#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/qmmachine.h>


namespace votca { namespace xtp {


class Log2Mps : public ctp::QMTool
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
    
    if(_package=="xtp"){
        throw std::runtime_error("XTP has no log file. For xtp package just run the partialcharges tool on you .orb file");
    }
    _logfile = opt->get(key+".logfile").as<string>();
    

    _mpsfile = (opt->exists(key+".mpsfile")) ? 
        opt->get(key+".mpsfile").as<string>() : "";
    if (_mpsfile == "") _mpsfile = _logfile.substr(0,_logfile.size()-4)+".mps";

    cout << endl << "... ... " << _logfile << " => " << _mpsfile << flush;
}


bool Log2Mps::Evaluate() {
    
    // Logger (required for QM package, so we can just as well use it)
    ctp::Logger log;
    log.setPreface(ctp::logINFO, "\n... ...");
    log.setPreface(ctp::logDEBUG, "\n... ...");
    log.setReportLevel(ctp::logDEBUG);
    log.setMultithreading(true);  
    
    // Set-up QM package
    
    CTP_LOG_SAVE(ctp::logINFO,log) << "Using package <" << _package << ">" << flush;
    QMPackage *qmpack = QMPackages().Create(_package);    
    qmpack->doGetCharges(true);
    qmpack->setLog(&log);
    qmpack->setRunDir(".");
    qmpack->setLogFileName(_logfile);
    
    // Create orbitals, fill with life & extract QM atoms
    Orbitals orbs;
    int cdx = qmpack->ParseLogFile(orbs);
    if (!cdx) {
        cout << "\nERROR Parsing " << _logfile << "failed. Abort." << endl;
        throw std::runtime_error("(see above, parsing error)");
    }    
    vector<QMAtom*> &qmatoms = orbs.QMAtoms();
    vector<QMAtom*>::iterator it;
    
    // Sanity checks, total charge
    double Q = 0.0;
    for (it = qmatoms.begin(); it < qmatoms.end(); ++it) {
        Q += (*it)->getPartialcharge();
    }
    
    if (qmatoms.size() < 1) {
        cout << "\nERROR No charges extracted from " << _logfile 
            << ". Abort.\n" << flush;
        throw std::runtime_error("(see above, input or parsing error)");
    }
    CTP_LOG_SAVE(ctp::logINFO,log) 
        << qmatoms.size() << " QM atoms, total charge Q = " << Q << flush;    
    
    
    // Convert to polar segment & write mps-file
    QMInterface qmmface;
    ctp::PolarSeg pseg = qmmface.Convert(qmatoms);
    
    string tag = "::LOG2MPS " 
        + (boost::format("(log-file='%1$s' : %2$d QM atoms)")
        % _logfile % qmatoms.size()).str();    
    pseg.WriteMPS(_mpsfile, tag);
    return true;
}



}}

#endif
