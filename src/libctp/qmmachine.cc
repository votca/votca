#include <votca/ctp/qmmachine.h>
#include <sys/stat.h>


namespace votca { namespace ctp {


QMMachine::QMMachine(XJob *job, XInductor *xind, Gaussian *qmpack,
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
    
    if (mav)
        printf("\n... ... ... dR %1.4f dQ %1.4f QM %1.4f MM %1.4f IT %d",
                _crit_dR, _crit_dQ, _crit_dE_QM, _crit_dE_MM, _maxIter);
}


QMMachine::~QMMachine() {
    
    vector<QMMIter*> ::iterator qit;
    for (qit = _iters.begin(); qit < _iters.end(); ++qit) {
        delete *qit;
    }
    _iters.clear();
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


bool QMMachine::Iterate(string jobFolder, int iterCnt) {

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

    printf("\n... ... ... Summary - iteration %d", iterCnt+1);
    printf("\n... ... ... ... QM Size = %d atoms", atoms.size());
    printf("\n... ... ... ... E(QM)   = %+4.9e", thisIter->getQMEnergy()); 
    printf("\n... ... ... ... E(SF)   = %+4.9e", thisIter->getSFEnergy());
    printf("\n... ... ... ... E(FM)   = %+4.9e", thisIter->getFMEnergy());
    printf("\n... ... ... ... E(MM)   = %+4.9e", thisIter->getMMEnergy());
    printf("\n... ... ... ... E(QMMM) = %+4.9e", thisIter->getQMMMEnergy());    
    printf("\n... ... ... ... RMS(dR) = %+4.9e", thisIter->getRMSdR());
    printf("\n... ... ... ... RMS(dQ) = %+4.9e", thisIter->getRMSdQ());
    printf("\n... ... ... ... SUM(dQ) = %+4.9e", thisIter->getSUMdQ());
    
    return 0;   
}


QMMachine::QMMIter *QMMachine::CreateNewIter() {
    
    QMMIter *newIter = new QMMIter(_iters.size());
    this->_iters.push_back(newIter);
    return newIter;
}


void QMMachine::WriteQMPackInputFile(string inputFile, Gaussian *qmpack, XJob *job) {
    
    // TODO _qmpack should do this entirely independently
    FILE *out;
    out = fopen(inputFile.c_str(), "w");
    _qmpack->WriteInputHeader(out, job->getTag());
    job->getPolarTop()->PrintInduState(out, _qmpack->getPackageName(), true, 1e-04);
    fclose(out);
    
}


bool QMMachine::hasConverged() {
    
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
    
    printf("\n... ... ... ... Convg dR = %s", _convg_dR ? "true" : "false");
    printf("\n... ... ... ... Convg dQ = %s", _convg_dQ ? "true" : "false");
    printf("\n... ... ... ... Convg QM = %s", _convg_dE_QM ? "true" : "false");
    printf("\n... ... ... ... Convg MM = %s", _convg_dE_MM ? "true" : "false");
    
    return _isConverged;
}


void QMMachine::QMMIter::ConvertPSitesToQMAtoms(vector< PolarSeg* > &psegs,
                                       vector< QMAtom * > &qmatoms) {
    
    assert(qmatoms.size() == 0);    
    return;   
}


void QMMachine::QMMIter::ConvertQMAtomsToPSites(vector< QMAtom* > &qmatoms,
                                       vector< PolarSeg* > &psegs) {
    assert(qmatoms.size() == 0);
    return;
}


void QMMachine::QMMIter::UpdatePosChrgFromQMAtoms(vector< QMAtom* > &qmatoms,
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


void QMMachine::QMMIter::setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM) {
    
    _hasdRdQ = true;    
    _dR_RMS = dR_RMS;
    _dQ_RMS = dQ_RMS;
    _dQ_SUM = dQ_SUM;
    return;
}


void QMMachine::QMMIter::setQMSF(double energy_QM, double energy_SF) {
    
    _hasQM = true;    
    _e_QM = energy_QM;
    _e_SF = energy_SF;    
    return;
}


void QMMachine::QMMIter::setE_FM(double ef00, double ef01, double ef02, 
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


double QMMachine::QMMIter::getMMEnergy() {
    
    assert(_hasMM);
    return _ef_11 + _ef_12 + _em_1_ + _em_2_;
}


double QMMachine::QMMIter::getQMMMEnergy() {
    
    assert(_hasQM && _hasMM);    
    return _e_QM + _ef_11 + _ef_12 + _em_1_ + _em_2_;    
}

   
    
    
    
    
    
}}
