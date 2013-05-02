#include <votca/ctp/qmmachine.h>
#include <sys/stat.h>
#include <boost/format.hpp>


using boost::format;

namespace votca { namespace ctp {

template<class QMPackage>
QMMachine<QMPackage>::QMMachine(XJob *job, XInductor *xind, QMPackage *qmpack,
                     Property *opt, string sfx, int nst, bool mav) 
                   : _job(job), _xind(xind), _qmpack(qmpack), _subthreads(nst),
                     _isConverged(false) {
    
    string key = sfx + ".qmmmconvg";
    if (opt->exists(key+".dR")) {
        _crit_dR = opt->get(key+".dR").as<double>();
    }
    else {
        _crit_dR = 0.01; // nm
    }
    if (opt->exists(key+".dQ")) {
        _crit_dQ = opt->get(key+".dQ").as<double>();
    }
    else {
        _crit_dQ = 0.01; // e
    }
    if (opt->exists(key+".dE_QM")) {
        _crit_dE_QM = opt->get(key+".dE_QM").as<double>();
    }
    else {
        _crit_dE_QM = 0.001; // eV
    }
    if (opt->exists(key+".dE_MM")) {
        _crit_dE_MM = opt->get(key+".dE_MM").as<double>();
    }
    else {
        _crit_dE_MM = _crit_dE_QM; // eV
    }
    if (opt->exists(key+".max_iter")) {
        _maxIter = opt->get(key+".max_iter").as<int>();
    }
    else {
        _maxIter = 32;
    }
}


template<class QMPackage>
QMMachine<QMPackage>::~QMMachine() {
    
    vector<QMMIter*> ::iterator qit;
    for (qit = _iters.begin(); qit < _iters.end(); ++qit) {
        delete *qit;
    }
    _iters.clear();
}


template<class QMPackage>
void QMMachine<QMPackage>::Evaluate(XJob *job) {
    
    LOG(logINFO,*_log)
       << format("... dR %1$1.4f dQ %2$1.4f QM %3$1.4f MM %4$1.4f IT %5$d")
       % _crit_dR % _crit_dQ % _crit_dE_QM % _crit_dE_MM % _maxIter << flush;
    
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
    LOG(logINFO,*_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;
    
    // SET ITERATION-TIME CONSTANTS
    _qmpack->setCharge(chrg);
    _qmpack->setSpin(spin);
    _qmpack->setThreads(_subthreads);
    
    int iterCnt = 0;
    int iterMax = _maxIter;
    for ( ; iterCnt < iterMax; ++iterCnt) {
        
        bool code = Iterate(jobFolder, iterCnt);
        if (hasConverged()) { break; }
    }
    
    if (iterCnt == iterMax-1 && !_isConverged) {
        printf("\n... ... ... WARNING: Not converged within %d iterations.", iterMax);
    }
    
    return;
}


template<class QMPackage>
bool QMMachine<QMPackage>::Iterate(string jobFolder, int iterCnt) {

    QMMIter *thisIter = this->CreateNewIter();
    int iter = iterCnt;
    string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);
    mkdir(runFolder.c_str(), 0755);    
    
    // RUN CLASSICAL INDUCTION & SAVE
    _xind->Evaluate(_job);
    assert(_xind->hasConverged());
    thisIter->setE_FM(_job->getEF00(), _job->getEF01(), _job->getEF02(),
                      _job->getEF11(), _job->getEF12(), _job->getEM0(),
                      _job->getEM1(),  _job->getEM2(),  _job->getETOT());
    
    // WRITE AND SET QM INPUT FILE
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
    
    // EXTRACT LOG-FILE INFOS TO ORBITALS
    Orbitals orb_iter;
    _qmpack->ParseLogFile(&orb_iter);
    
    assert(orb_iter.hasSelfEnergy());
    assert(orb_iter.hasQMEnergy());
    
    // EXTRACT & SAVE QM ENERGIES
    double energy___sf = orb_iter.getSelfEnergy();
    double energy_qmsf = orb_iter.getQMEnergy();
    double energy_qm__ = energy_qmsf - energy___sf;
    thisIter->setQMSF(energy_qm__, energy___sf);
    _job->setEnergy_QMMM(thisIter->getQMEnergy(), thisIter->getSFEnergy(),
                         thisIter->getQMMMEnergy());
    
    // EXTRACT & SAVE QMATOM DATA
    vector< QMAtom* > &atoms = *(orb_iter.getAtoms());    
    thisIter->UpdatePosChrgFromQMAtoms(atoms, _job->getPolarTop()->QM0());

    LOG(logINFO,*_log) 
        << format("Summary - iteration %1$d:") % (iterCnt+1) << flush;
    LOG(logINFO,*_log)
        << format("... QM Size  = %1$d atoms") % int(atoms.size()) << flush;
    LOG(logINFO,*_log)
        << format("... E(QM)    = %1$+4.9e") % thisIter->getQMEnergy() << flush;
    LOG(logINFO,*_log)
        << format("... E(SF)    = %1$+4.9e") % thisIter->getSFEnergy() << flush;
    LOG(logINFO,*_log)
        << format("... E(FM)    = %1$+4.9e") % thisIter->getFMEnergy() << flush;
    LOG(logINFO,*_log)
        << format("... E(MM)    = %1$+4.9e") % thisIter->getMMEnergy() << flush;
    LOG(logINFO,*_log)
        << format("... E(QMMM)  = %1$+4.9e") % thisIter->getQMMMEnergy() << flush;
    LOG(logINFO,*_log)
        << format("... RMS(dR)  = %1$+4.9e") % thisIter->getRMSdR() << flush;
    LOG(logINFO,*_log)
        << format("... MS(dQ)   = %1$+4.9e") % thisIter->getRMSdQ() << flush;
    LOG(logINFO,*_log)
        << format("... SUM(dQ)  = %1$+4.9e") % thisIter->getSUMdQ() << flush;
    
    return 0;
}


template<class QMPackage>
QMMIter *QMMachine<QMPackage>::CreateNewIter() {
    
    QMMIter *newIter = new QMMIter(_iters.size());
    this->_iters.push_back(newIter);
    return newIter;
}


template<class QMPackage>
void QMMachine<QMPackage>::WriteQMPackInputFile(string inputFile, QMPackage *qmpack, XJob *job) {
    
    // TODO _qmpack should do this entirely independently
    FILE *out;
    out = fopen(inputFile.c_str(), "w");
    _qmpack->WriteInputHeader(out, job->getTag());
    job->getPolarTop()->PrintInduState(out, _qmpack->getPackageName(), true, 1e-04);
    fclose(out);
    
}


template<class QMPackage>
bool QMMachine<QMPackage>::hasConverged() {
    
    _convg_dR = false;
    _convg_dQ = false;
    _convg_dE_QM = false;
    _convg_dE_MM = false;
    
    if (_iters.size() > 1) {
        
        QMMIter *iter_0 = _iters[_iters.size()-2];
        QMMIter *iter_1 = _iters[_iters.size()-1];
        
        double dR = iter_1->getRMSdR();
        double dQ = iter_1->getRMSdQ();
        double dE_QM = iter_1->getQMEnergy() - iter_0->getQMEnergy();
        double dE_MM = iter_1->getMMEnergy() - iter_0->getMMEnergy();
        
        if (dR <= _crit_dR) _convg_dR = true;
        if (dQ <= _crit_dQ) _convg_dQ = true;
        if (dE_QM*dE_QM <= _crit_dE_QM*_crit_dE_QM) _convg_dE_QM = true;
        if (dE_MM*dE_MM <= _crit_dE_MM*_crit_dE_MM) _convg_dE_MM = true;        
    }
    
    _isConverged = ((_convg_dR && _convg_dQ) && (_convg_dE_QM && _convg_dE_MM));
    
    LOG(logINFO,*_log) 
        << format("... Convg dR = %s") % (_convg_dR ? "true" : "false") << flush;
    LOG(logINFO,*_log) 
        << format("... Convg dQ = %s") % (_convg_dQ ? "true" : "false") << flush;
    LOG(logINFO,*_log) 
        << format("... Convg QM = %s") % (_convg_dE_QM ? "true" : "false") << flush;
    LOG(logINFO,*_log) 
        << format("... Convg MM = %s") % (_convg_dE_MM ? "true" : "false") << flush;
    
    return _isConverged;
}


void QMMIter::ConvertPSitesToQMAtoms(vector< PolarSeg* > &psegs,
                                       vector< QMAtom * > &qmatoms) {
    
    assert(qmatoms.size() == 0);    
    return;   
}


void QMMIter::ConvertQMAtomsToPSites(vector< QMAtom* > &qmatoms,
                                       vector< PolarSeg* > &psegs) {
    assert(qmatoms.size() == 0);
    return;
}


void QMMIter::UpdatePosChrgFromQMAtoms(vector< QMAtom* > &qmatoms,
                                         vector< PolarSeg* > &psegs) {
    
    double AA_to_NM = 0.1; // Angstrom to nanometer
    
    double dR_RMS = 0.0;
    double dQ_RMS = 0.0;
    double dQ_SUM = 0.0;
    
    for (int i = 0, qac = 0; i < psegs.size(); ++i) {
        PolarSeg *pseg = psegs[i];
        for (int j = 0; j < pseg->size(); ++j, ++qac) {
            
            // Retrieve info from QMAtom
            QMAtom *qmatm = qmatoms[qac];
            vec upd_r = vec(qmatm->x, qmatm->y, qmatm->z);
            upd_r *= AA_to_NM;
            double upd_Q00 = qmatm->charge;
            
            // Compare to previous r, Q00
            APolarSite *aps = (*pseg)[j];
            vec old_r = aps->getPos();
            double old_Q00 = aps->getQ00();
            double dR = abs(upd_r - old_r);
            double dQ00 = upd_Q00 - old_Q00;
            
            dR_RMS += dR*dR;
            dQ_RMS += dQ00*dQ00;
            dQ_SUM += dQ00;
            
            // Forward updated r, Q00 to APS
            aps->setPos(upd_r);
            aps->setQ00(upd_Q00, 0);            
        }
    }
    
    dR_RMS /= qmatoms.size();
    dQ_RMS /= qmatoms.size();
    dR_RMS = sqrt(dR_RMS);
    dQ_RMS = sqrt(dQ_RMS);

    this->setdRdQ(dR_RMS, dQ_RMS, dQ_SUM);
}


void QMMIter::setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM) {
    
    _hasdRdQ = true;    
    _dR_RMS = dR_RMS;
    _dQ_RMS = dQ_RMS;
    _dQ_SUM = dQ_SUM;
    return;
}


void QMMIter::setQMSF(double energy_QM, double energy_SF) {
    
    _hasQM = true;    
    _e_QM = energy_QM;
    _e_SF = energy_SF;    
    return;
}


void QMMIter::setE_FM(double ef00, double ef01, double ef02, 
    double ef11, double ef12, double em0, double em1,  double em2, double efm) {
    
    _hasMM = true;
    _ef_00 = ef00;
    _ef_01 = ef01;
    _ef_02 = ef02;
    _ef_11 = ef11;
    _ef_12 = ef12;
    _em_0_ = em0;
    _em_1_ = em1;
    _em_2_ = em2;
    _e_fm_ = efm;
    return;
}


double QMMIter::getMMEnergy() {
    
    assert(_hasMM);
    return _ef_11 + _ef_12 + _em_1_ + _em_2_;
}


double QMMIter::getQMMMEnergy() {
    
    assert(_hasQM && _hasMM);    
    return _e_QM + _ef_11 + _ef_12 + _em_1_ + _em_2_;    
}


// REGISTER QM PACKAGES
template class QMMachine<Gaussian>;
    
    
    
}}
