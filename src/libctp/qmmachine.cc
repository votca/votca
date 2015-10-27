// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/qmmachine.h>
#include <sys/stat.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/ctp/logger.h>
#include <votca/ctp/elements.h>
#include <votca/tools/linalg.h>
#include <votca/ctp/espfit.h>

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
    
    key = sfx + ".control";
    if (opt->exists(key+".split_dpl")) {
        _split_dpl = opt->get(key+".split_dpl").as<bool>();
    }
    else {
        _split_dpl = true;
    }
    if (opt->exists(key+".dpl_spacing")) {
        _dpl_spacing = opt->get(key+".dpl_spacing").as<double>();
    }
    else {
        _dpl_spacing = 1e-3;
    }
    
    // GWBSE options
    key = sfx + ".gwbse";
    string _gwbse_xml = opt->get(key + ".gwbse_options").as<string> ();
    //cout << endl << "... ... Parsing " << _package_xml << endl ;
    load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());
    _state = opt->get(key+".state").as< int >();
    _type  = opt->get(key+".type").as< string >();
    if ( _type != "singlet" && _type != "triplet" ){
        throw runtime_error(" Invalid excited state type! " + _type);
    }
   
    key = sfx + ".gwbse.filter";
    if (opt->exists(key + ".oscillator_strength") && _type != "triplet" ){
        _has_osc_filter = true;
        _osc_threshold  =  opt->get(key + ".oscillator_strength").as<double> ();
    } else {
        _has_osc_filter = false;
    }
        
    if (opt->exists(key + ".charge_transfer")  ){
        _has_dQ_filter = true;
        _dQ_threshold  =  opt->get(key + ".charge_transfer").as<double> ();
    } else {
        _has_dQ_filter = false;
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
    
    // FIGURE OUT CHARGE + MULTIPLICITY
    double dQ = 0.0;
    for (unsigned int i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
        dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
    }
    int chrg = round(dQ);
    int spin = ( (chrg < 0) ? -chrg:chrg ) % 2 + 1;
    LOG(logINFO,*_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;
    
    
    // PREPARE JOB DIRECTORY
    string jobFolder = "xjob_" + boost::lexical_cast<string>(_job->getId())
                     + "_" + _job->getTag();    
    bool created = boost::filesystem::create_directory(jobFolder);
    if (created){ 
        LOG(logINFO,*_log) << "Created directory " << jobFolder << flush;
    }
    
    
    // SET ITERATION-TIME CONSTANTS
    // TO ADJUST
    
    _qmpack->setCharge(chrg);
    _qmpack->setSpin(spin);

    //_qmpack->setThreads(_subthreads);
    

    int iterCnt = 0;
    int iterMax = _maxIter;
    for ( ; iterCnt < iterMax; ++iterCnt) {
        
        //bool code = 
	(void)Iterate(jobFolder, iterCnt);
        if (hasConverged()) { break; }
    }
    
    if (iterCnt == iterMax-1 && !_isConverged) {
        LOG(logWARNING,*_log)
            << format("Not converged within %1$d iterations.") % iterMax;
    }
    
    return;
}


template<class QMPackage>
bool QMMachine<QMPackage>::Iterate(string jobFolder, int iterCnt) {

    // CREATE ITERATION OBJECT & SETUP RUN DIRECTORY
    QMMIter *thisIter = this->CreateNewIter();
    int iter = iterCnt;
    string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);
       
    bool created = boost::filesystem::create_directory(runFolder);
    if (created) 
        LOG(logDEBUG,*_log) << "Created directory " << runFolder << flush;
    else
        LOG(logWARNING,*_log) << "Could not create directory " << runFolder << flush;
    
    
    // RUN CLASSICAL INDUCTION & SAVE
    _job->getPolarTop()->PrintPDB(runFolder + "/QM0_MM1_MM2.pdb");
    _xind->Evaluate(_job);
    assert(_xind->hasConverged());
    thisIter->setE_FM(_job->getEF00(), _job->getEF01(), _job->getEF02(),
                      _job->getEF11(), _job->getEF12(), _job->getEM0(),
                      _job->getEM1(),  _job->getEM2(),  _job->getETOT());
    
    // WRITE AND SET QM INPUT FILE
    Orbitals orb_iter_input;
    
    vector<Segment*> empty;
    thisIter->GenerateQMAtomsFromPolarSegs(_job->getPolarTop(), orb_iter_input, _split_dpl, _dpl_spacing);
      
    _qmpack->setRunDir(runFolder);
    
    LOG(logDEBUG,*_log) << "Writing input file " << runFolder << flush;
    
    _qmpack->WriteInputFile(empty, &orb_iter_input);
 
    FILE *out;
    out = fopen((runFolder + "/system.pdb").c_str(),"w");
    orb_iter_input.WritePDB( out );
    fclose(out);
         
    // RUN HERE (OVERRIDE - COPY EXISTING LOG-FILE)
    //string cpstr = "cp e_1_n.log " + path_logFile;
    //int sig = system(cpstr.c_str());
    //_qmpack->setLogFileName(path_logFile);
    
    //Commented out for test Jens 
    _qmpack->Run();
    
    // EXTRACT LOG-FILE INFOS TO ORBITALS   
    Orbitals orb_iter_output;
    _qmpack->ParseLogFile(&orb_iter_output);
    
    
    
     // Ground state density matrix
    // ub::matrix<double> &_dft_orbitals_GS = orb_iter_output.MOCoefficients();
    // int _parse_orbitals_status_GS = _qmpack->ParseOrbitalsFile( &orb_iter_output );


    
    
    /* // AOESP matrix test
    // load DFT basis set (element-wise information) from xml file
    BasisSet dftbs;
    //dftbs.LoadBasisSet( orb_iter_output.getDFTbasis() );
    dftbs.LoadBasisSet( "ubecppol" );
    //LOG(logDEBUG, *_log) << TimeStamp() << " Loaded DFT Basis Set " <<  orb_iter_output.getDFTbasis()  << flush;
    
    // fill DFT AO basis by going through all atoms 
    AOBasis dftbasis;
    dftbasis.AOBasisFill(&dftbs, orb_iter_output.QMAtoms() );
    dftbasis.ReorderMOs(_dft_orbitals_GS, orb_iter_output.getQMpackage(), "votca" );
    ub::matrix<double> &DMATGS=orb_iter_output.DensityMatrixGroundState(_dft_orbitals_GS);
    vector< QMAtom* >& Atomlist= orb_iter_output.QMAtoms();
    
    Espfit esp;
    // Espfit esp
    //Espfit esp(vector< QMAtom* >& Atomlist, ub::matrix<double> &DMATGS, AOBasis &dftbasis);
    esp.setLog(_log);
    esp.FittoDensity(Atomlist, DMATGS, dftbasis); */
  
    

 
    
    

    // GW-BSE starts here
    bool _do_gwbse = true; // needs to be set by options!!!
    double energy___ex = 0.0;
    
    if ( _do_gwbse ){


        // for GW-BSE, we also need to parse the orbitals file
        _qmpack->ParseOrbitalsFile( &orb_iter_output );
        //int _parse_orbitals_status = _qmpack->ParseOrbitalsFile( &orb_iter_output );
        std::vector<int> _state_index;
       _gwbse.Initialize( &_gwbse_options );
      if ( _state > 0 ){
        LOG(logDEBUG,*_log) << "Excited state via GWBSE: " <<  flush;
        LOG(logDEBUG,*_log) << "  --- type:              " << _type << flush;
        LOG(logDEBUG,*_log) << "  --- state:             " << _state << flush;
        if ( _has_osc_filter) { LOG(logDEBUG,*_log) << "  --- filter: osc.str. > " << _osc_threshold << flush;}
        if ( _has_dQ_filter)  {LOG(logDEBUG,*_log) << "  --- filter: crg.trs. > " << _dQ_threshold << flush;}
        
        if ( _has_osc_filter && _has_dQ_filter ){
            LOG(logDEBUG,*_log) << "  --- WARNING: filtering for optically active CT transition - might not make sense... "  << flush;
        }
        
        // define own logger for GW-BSE that is written into a runFolder logfile
        Logger gwbse_logger(logDEBUG);
        gwbse_logger.setMultithreading(false);
        _gwbse.setLogger(&gwbse_logger);
        gwbse_logger.setPreface(logINFO,    (format("\nGWBSE INF ...") ).str());
        gwbse_logger.setPreface(logERROR,   (format("\nGWBSE ERR ...") ).str());
        gwbse_logger.setPreface(logWARNING, (format("\nGWBSE WAR ...") ).str());
        gwbse_logger.setPreface(logDEBUG,   (format("\nGWBSE DBG ...") ).str());
        
        // actual GW-BSE run
        _gwbse.Evaluate( &orb_iter_output );
        //bool _evaluate = _gwbse.Evaluate( &orb_iter_output );
        
       
        // write logger to log file
        ofstream ofs;
        string gwbse_logfile = runFolder + "/gwbse.log";
        ofs.open(gwbse_logfile.c_str(), ofstream::out);
        if (!ofs.is_open()) {
            throw runtime_error("Bad file handle: " + gwbse_logfile);
        }    
        ofs << gwbse_logger << endl;
        ofs.close();

        // PROCESSING the GW-BSE result
        // - find the excited state of interest
        // oscillator strength filter

        
        if ( _has_osc_filter ){
            
            // go through list of singlets
            const std::vector<std::vector<double> >& TDipoles = orb_iter_output.TransitionDipoles();
            for (unsigned _i=0; _i < TDipoles.size(); _i++ ) {
                
                double osc = (TDipoles[_i][0] * TDipoles[_i][0] + TDipoles[_i][1] * TDipoles[_i][1] + TDipoles[_i][2] * TDipoles[_i][2]) * 1.0 / 3.0 * (orb_iter_output.BSESingletEnergies()[_i]) ;
                if ( osc > _osc_threshold ) _state_index.push_back(_i);
            } 
            
      
            
        } else {
            
            if ( _type == "singlet" ){
               for (unsigned _i=0; _i < orb_iter_output.TransitionDipoles().size(); _i++ ) {
                   _state_index.push_back(_i);
               }
            } else {
               for (unsigned _i=0; _i < orb_iter_output.BSETripletEnergies().size(); _i++ ) {
                   _state_index.push_back(_i);
               }
            }
        }


        // filter according to charge transfer, go through list of excitations in _state_index
         if  (_has_dQ_filter ) {
            std::vector<int> _state_index_copy;
            if ( _type == "singlets" ){
             // go through list of singlets
            const std::vector<double>& dQ_fragA = orb_iter_output.FragmentAChargesSingEXC();
            //const std::vector<double>& dQ_fragB = orb_iter_output.FragmentBChargesSingEXC();
            for (unsigned _i=0; _i < _state_index.size(); _i++ ) {
                if ( std::abs(dQ_fragA[_i]) > _dQ_threshold ) {
                    _state_index_copy.push_back(_state_index[_i]);
                }
            } 
            _state_index = _state_index_copy;
            } else if ( _type == "triplets"){
              // go through list of triplets
            const std::vector<double>& dQ_fragA = orb_iter_output.FragmentAChargesTripEXC();
            //const std::vector<double>& dQ_fragB = orb_iter_output.FragmentBChargesTripEXC();
            for (unsigned _i=0; _i < _state_index.size(); _i++ ) {
                if ( std::abs(dQ_fragA[_i]) > _dQ_threshold ) {
                    _state_index_copy.push_back(_state_index[_i]);
                }
            } 
            _state_index = _state_index_copy;               
                
                
            }
         }
        
        
        if ( _state_index.size() < 1 ){
            throw runtime_error("Excited state filter yields no states! ");
            
        }
        // - output its energy
        if ( _type == "singlet" ){
            energy___ex = orb_iter_output.BSESingletEnergies()[_state_index[_state-1]]*13.6058; // to eV
        } else if ( _type == "triplet" ) {
            energy___ex = orb_iter_output.BSETripletEnergies()[_state_index[_state-1]]*13.6058; // to eV
        }
   
        // ub::matrix<double> &_dft_orbitals_GS = orb_iter_output.MOCoefficients();
        // int _parse_orbitals_status_GS = _qmpack->ParseOrbitalsFile( &orb_iter_output );

      } // only if state >0
        
        // calculate density matrix for this excited state
        ub::matrix<double> &_dft_orbitals = orb_iter_output.MOCoefficients();
        // load DFT basis set (element-wise information) from xml file
        BasisSet dftbs;
        if ( orb_iter_output.getDFTbasis() != "" ) {
          dftbs.LoadBasisSet( orb_iter_output.getDFTbasis() );
	} else{
	  dftbs.LoadBasisSet( _gwbse.get_dftbasis_name() );

	}  
       LOG(logDEBUG, *_log) << TimeStamp() << " Loaded DFT Basis Set " <<  orb_iter_output.getDFTbasis()  << flush;

    
    
    
    
        // fill DFT AO basis by going through all atoms 
        AOBasis dftbasis;
        dftbasis.AOBasisFill(&dftbs, orb_iter_output.QMAtoms() );
        dftbasis.ReorderMOs(_dft_orbitals, orb_iter_output.getQMpackage(), "votca" );
        // TBD: Need to switch between singlets and triplets depending on _type
        ub::matrix<double> &DMATGS=orb_iter_output.DensityMatrixGroundState(_dft_orbitals);

        ub::matrix<double> DMAT_tot=DMATGS; // Ground state + hole_contribution + electron contribution

	if ( _state > 0 ){ 
	  ub::matrix<float>& BSECoefs = orb_iter_output.BSESingletCoefficients();
	  std::vector<ub::matrix<double> > &DMAT = orb_iter_output.DensityMatrixExcitedState( _dft_orbitals , BSECoefs, _state_index[_state-1]);
	  DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1]; // Ground state + hole_contribution + electron contribution
	}
    
        // fill DFT AO basis by going through all atoms 
        vector< QMAtom* >& Atomlist= orb_iter_output.QMAtoms();
        

        
        Espfit esp;
    
        esp.setLog(_log);
        esp.Fit2Density(Atomlist, DMAT_tot, dftbasis,dftbs,"medium",false);
    
    
    
    
       

    }
    
    

    out = fopen((runFolder + "/parsed.pdb").c_str(),"w");
    orb_iter_input.WritePDB( out );
    fclose(out);
    
    assert(orb_iter_output.hasSelfEnergy());
    assert(orb_iter_output.hasQMEnergy());
    
    // EXTRACT & SAVE QM ENERGIES
    double energy___sf = orb_iter_output.getSelfEnergy();
    double energy_qmsf = orb_iter_output.getQMEnergy();
    double energy_qm__ = energy_qmsf - energy___sf ;
    thisIter->setQMSF(energy_qm__, energy___sf, energy___ex);
    _job->setEnergy_QMMM(thisIter->getQMEnergy(),thisIter->getGWBSEEnergy(), thisIter->getSFEnergy(),
                         thisIter->getQMMMEnergy());
    
    // EXTRACT & SAVE QMATOM DATA
    vector< QMAtom* > &atoms = *(orb_iter_output.getAtoms());
    
    thisIter->UpdatePosChrgFromQMAtoms(atoms, _job->getPolarTop()->QM0());

    LOG(logINFO,*_log) 
        << format("Summary - iteration %1$d:") % (iterCnt+1) << flush;
    LOG(logINFO,*_log)
        << format("... QM Size  = %1$d atoms") % int(atoms.size()) << flush;
    LOG(logINFO,*_log)
        << format("... E(QM)    = %1$+4.9e") % thisIter->getQMEnergy() << flush;
    LOG(logINFO,*_log)
        << format("... E(GWBSE) = %1$+4.9e") % thisIter->getGWBSEEnergy() << flush;
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
        << format("... RMS(dQ)  = %1$+4.9e") % thisIter->getRMSdQ() << flush;
    LOG(logINFO,*_log)
        << format("... SUM(dQ)  = %1$+4.9e") % thisIter->getSUMdQ() << flush;
    
    // CLEAN DIRECTORY
    _qmpack->CleanUp();

    
    /*
    int removed = boost::filesystem::remove_all(runFolder);
    if (removed > 0) 
        LOG(logDEBUG,*_log) << "Removed directory " << runFolder << flush;
    else 
        LOG(logWARNING,*_log) << "Could not remove dir " << runFolder << flush;
    */
    return 0;
     
}


template<class QMPackage>
QMMIter *QMMachine<QMPackage>::CreateNewIter() {
    
    QMMIter *newIter = new QMMIter(_iters.size());
    this->_iters.push_back(newIter);
    return newIter;
}

/*
template<class QMPackage>
void QMMachine<QMPackage>::WriteQMPackInputFile(string inputFile, QMPackage *qmpack, XJob *job) {
    
    // TODO _qmpack should do this entirely independently
    FILE *out;
    out = fopen(inputFile.c_str(), "w");

    // TO ADJUST
    //_qmpack->WriteInputHeader(out, job->getTag());
    job->getPolarTop()->PrintInduState(out, _qmpack->getPackageName(), true, 1e-04);
    fclose(out);
    
}
*/

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
    
    for (unsigned int i = 0, qac = 0; i < psegs.size(); ++i) {
        PolarSeg *pseg = psegs[i];
        for (unsigned int j = 0; j < pseg->size(); ++j, ++qac) {
            
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


void QMMIter::GenerateQMAtomsFromPolarSegs(PolarTop *ptop, Orbitals &orb, 
        bool split_dpl, double dpl_spacing) {
    
    double AA_to_NM = 0.1; // Angstrom to nanometer
    
    // INNER SHELL QM0
    for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
        PolarSeg *pseg = ptop->QM0()[i];
        for (unsigned int j = 0; j < pseg->size(); ++j) {
            
            APolarSite *aps = (*pseg)[j];
            vec pos = aps->getPos()/AA_to_NM;
            double Q = aps->getQ00();
            string type = "qm";

            orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, false);            
              
        }
    }
    
    // MIDDLE SHELL MM1
    for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
        PolarSeg *pseg = ptop->MM1()[i];
        for (unsigned int j = 0; j < pseg->size(); ++j) {
            
            APolarSite *aps = (*pseg)[j];
            vec pos = aps->getPos()/AA_to_NM;
            double Q = aps->getQ00();
            string type = "mm";

            orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, true);
            
            if (split_dpl) {
                vec tot_dpl = vec(aps->U1x,aps->U1y,aps->U1z);
                if (aps->getRank() > 0)
                    { tot_dpl += vec(aps->Q1x,aps->Q1y,aps->Q1z); }            
                // Calculate virtual charge positions
                double a        = dpl_spacing; // this is in nm
                double mag_d    = abs(tot_dpl); // this is in e * nm
                vec    dir_d_0  = tot_dpl.normalize(); 
                vec    dir_d    = dir_d_0.normalize();
                vec    A        = pos + 0.5 * a * dir_d /AA_to_NM; // converted to AA
                vec    B        = pos - 0.5 * a * dir_d /AA_to_NM;
                double qA       = mag_d / a;
                double qB       = - qA;
                // Zero out if magnitude small [e*nm]
                if (aps->getIsoP() < 1e-9 || mag_d < 1e-9) {
                    A = aps->getPos() + 0.1*a*vec(1,0,0); // != pos since self-energy may diverge
                    B = aps->getPos() - 0.1*a*vec(1,0,0);
                    qA = 0;
                    qB = 0;
                }
                orb.AddAtom("A", A.x(), A.y(), A.z(), qA, true);
                orb.AddAtom("B", B.x(), B.y(), B.z(), qB, true);
            }             
        }
    }
    
    // OUTER SHELL MM2
    for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
        PolarSeg *pseg = ptop->MM2()[i];
        for (unsigned int j = 0; j < pseg->size(); ++j) {
            
            APolarSite *aps = (*pseg)[j];
            vec pos = aps->getPos()/AA_to_NM;
            double Q = aps->getQ00();
            string type = "mm";

            orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, true);              
        }
    }
    return;
    
    
}


void QMMIter::setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM) {
    
    _hasdRdQ = true;    
    _dR_RMS = dR_RMS;
    _dQ_RMS = dQ_RMS;
    _dQ_SUM = dQ_SUM;
    return;
}


void QMMIter::setQMSF(double energy_QM, double energy_SF, double energy_GWBSE) {
    
    _hasQM = true;
    _e_QM = energy_QM;
    _e_SF = energy_SF;    

    _hasGWBSE = true;
    _e_GWBSE = energy_GWBSE;
   
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
    
    assert(_hasQM && _hasMM && _hasGWBSE);    
    return _e_QM + + _e_GWBSE + _ef_11 + _ef_12 + _em_1_ + _em_2_;    
}


// REGISTER QM PACKAGES
template class QMMachine<QMPackage>;
    
    
    
}}
