/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/xtp/qmapemachine.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/ctp/logger.h>
#include <votca/xtp/elements.h>
#include <votca/xtp/espfit.h>

using boost::format;

namespace votca { namespace xtp {

template<class QMPackage>
QMAPEMachine<QMPackage>::QMAPEMachine(ctp::XJob *job, ctp::Ewald3DnD *cape, QMPackage *qmpack,
	 Property *opt, string sfx, int nst)
   : _subthreads(nst),_job(job), _qmpack(qmpack), _cape(cape),_isConverged(false) {
    
	// CONVERGENCE THRESHOLDS
    string key = sfx + ".convergence";
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

	// TASKS
	key = sfx + ".tasks";
		_run_ape = opt->get(key+".run_ape").as<bool>();
		_run_dft = opt->get(key+".run_dft").as<bool>();
		_run_gwbse = opt->get(key+".run_gwbse").as<bool>();
        
        
    // FITTING GRIDS
   // key = sfx+ ".grids";    
        


	// GWBSE CONFIG
    key = sfx + ".gwbse";
		string gwbse_xml = opt->get(key + ".gwbse_options").as<string>();
		load_property_from_xml(_gwbse_options, gwbse_xml.c_str());

		_state = opt->get(key+".state").as<int>();
		_type  = opt->get(key+".type").as<string>();
		if ( _type != "singlet" && _type != "triplet") {
			throw runtime_error(" Invalid excited state type! " + _type);
		}

		key = sfx + ".gwbse.filter";
		if (opt->exists(key + ".oscillator_strength") && _type != "triplet") {
			_has_osc_filter = true;
			_osc_threshold  =  opt->get(key + ".oscillator_strength").as<double>();
		}
		else {
			_has_osc_filter = false;
		}
		if (opt->exists(key + ".charge_transfer")) {
			_has_dQ_filter = true;
			_dQ_threshold  =  opt->get(key + ".charge_transfer").as<double>();
		}
		else {
			_has_dQ_filter = false;
		}
}


template<class QMPackage>
QMAPEMachine<QMPackage>::~QMAPEMachine() {
    
    std::vector<QMMIter*> ::iterator qit;
    for (qit = _iters.begin(); qit < _iters.end(); ++qit) {
        delete *qit;
    }
    _iters.clear();
}


template<class QMPackage>
void QMAPEMachine<QMPackage>::Evaluate(ctp::XJob *job) {
    
	// PREPARE JOB DIRECTORY
	string jobFolder = "job_" + boost::lexical_cast<string>(_job->getId())
					 + "_" + _job->getTag();
	bool created = boost::filesystem::create_directory(jobFolder);

	LOG(ctp::logINFO,*_log) << flush;
	if (created) {
		LOG(ctp::logINFO,*_log) << "Created directory " << jobFolder << flush;
        }

    LOG(ctp::logINFO,*_log)
       << format("... dR %1$1.4f dQ %2$1.4f QM %3$1.4f MM %4$1.4f IT %5$d")
       % _crit_dR % _crit_dQ % _crit_dE_QM % _crit_dE_MM % _maxIter << flush;
    
    // FIGURE OUT CHARGE + MULTIPLICITY
    double dQ = 0.0;
    for (unsigned i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
        dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
    }
    int chrg = round(dQ);
    int spin = ( (chrg < 0) ? -chrg:chrg ) % 2 + 1;
    LOG(ctp::logINFO,*_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;

    // SET ITERATION-TIME CONSTANTS
    _qmpack->setCharge(chrg);
    _qmpack->setSpin(spin);

  
    
 
    
  
    
    int iterCnt = 0;
    int iterMax = _maxIter;
    for ( ; iterCnt < iterMax; ++iterCnt) {
        
        //bool code = Iterate(jobFolder, iterCnt);
        Iterate(jobFolder, iterCnt);
        if (hasConverged()) { break; }
    }
    
    if (iterCnt == iterMax-1 && !_isConverged) {
        LOG(ctp::logWARNING,*_log)
            << format("Not converged within %1$d iterations.") % iterMax;
    }
    
    return;
}


template<class QMPackage>
bool QMAPEMachine<QMPackage>::Iterate(string jobFolder, int iterCnt) {

    // CREATE ITERATION OBJECT & SETUP RUN DIRECTORY
    QMMIter *thisIter = this->CreateNewIter();
    int iter = iterCnt;
    string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);
       
    LOG(ctp::logINFO,*_log) << flush;
    bool created = boost::filesystem::create_directory(runFolder);
    if (created) 
        LOG(ctp::logDEBUG,*_log) << "Created directory " << runFolder << flush;
    else
        LOG(ctp::logWARNING,*_log) << "Could not create directory " << runFolder << flush;

    // COMPUTE POLARIZATION STATE WITH QM0(0)
    if (_run_ape) {
		if (iterCnt == 0) {
			_cape->ShowAgenda(_log);
			// Reset FGC, start from BGP state, apply FP fields (BG & FG)
			_cape->EvaluateInductionQMMM(true, true, true, true, true);
		}
    
    
        
		std::vector< ctp::PolarSeg* > target_bg;     
                std::vector< ctp::PolarSeg* > target_fg;     
        
       
		if (iterCnt == 0) {
			// Add BG, do not add MM1 & QM0
			_cape->EvaluatePotential(target_bg, true, false, false);
		}
		// Do not add BG & QM0, add MM1
		_cape->EvaluatePotential(target_fg, false, true, false);
    }
    
    

    // COMPUTE WAVEFUNCTION & QM ENERGY
    // Generate charge shell from potentials
    
    
    
   
       
    
    
    // Run DFT
    Orbitals orb_iter_input;
    std::vector<ctp::Segment*> empty;
   
   
	_qmpack->setRunDir(runFolder);
	_qmpack->WriteInputFile(empty, &orb_iter_input);

	FILE *out;
	out = fopen((runFolder + "/system.pdb").c_str(),"w");
	orb_iter_input.WritePDB(out);
	fclose(out);

	Orbitals orb_iter_output;
	if (_run_dft) {
		_qmpack->Run();
	}
	_qmpack->ParseLogFile(&orb_iter_output);


    // Run GWBSE
	if (_run_gwbse){
		this->EvaluateGWBSE(orb_iter_output, runFolder);
	}

	// COMPUTE POLARIZATION STATE WITH QM0(n+1)
	if (_run_ape) {
		// Update QM0 density: QM0(n) => QM0(n+1)
		// ...
        thisIter->UpdatePosChrgFromQMAtoms(orb_iter_output.QMAtoms(),
            _job->getPolarTop()->QM0());

		// Do not reset FGC (= only reset FU), do not use BGP state, nor apply FP fields (BG & FG)
		_cape->EvaluateInductionQMMM(false, false, false, false, false);

		// COMPUTE MM ENERGY
		_cape->EvaluateEnergyQMMM();
		_cape->ShowEnergySplitting(_log);
	}

	// COMPILE HAMILTONIAN & CHECK FOR CONVERGENCE
	// ...


    // THIS IS A SHORT VERSION OF PEWALD3D WITH QM/MM ENERGY SPLITTING
    _cape->ShowAgenda(_log);
    _cape->EvaluateInductionQMMM(true, true, true, true, true);
    _cape->EvaluateEnergyQMMM();
    _cape->ShowEnergySplitting(_log);

    return true;



    // RUN CLASSICAL INDUCTION & SAVE
    _job->getPolarTop()->PrintPDB(runFolder + "/QM0_MM1_MM2.pdb");
    _xind->Evaluate(_job);
    assert(_xind->hasConverged());
    thisIter->setE_FM(_job->getEF00(), _job->getEF01(), _job->getEF02(),
                      _job->getEF11(), _job->getEF12(), _job->getEM0(),
                      _job->getEM1(),  _job->getEM2(),  _job->getETOT());
    
    // WRITE AND SET QM INPUT FILE

    _qmpack->setRunDir(runFolder);
    
    LOG(ctp::logDEBUG,*_log) << "Writing input file " << runFolder << flush;
    
    _qmpack->WriteInputFile(empty, &orb_iter_input);
         
    // RUN HERE (OVERRIDE - COPY EXISTING LOG-FILE)
    //string cpstr = "cp e_1_n.log " + path_logFile;
    //int sig = std::system(cpstr.c_str());
    //_qmpack->setLogFileName(path_logFile);
    
    //Commented out for test Jens 
    _qmpack->Run();
    
    // EXTRACT LOG-FILE INFOS TO ORBITALS   

    _qmpack->ParseLogFile(&orb_iter_output);
    
   
    // GW-BSE starts here
    bool _do_gwbse = true; // needs to be set by options!!!
    
    if (_do_gwbse){
    	this->EvaluateGWBSE(orb_iter_output, runFolder);
    }
    
    

    out = fopen((runFolder + "/parsed.pdb").c_str(),"w");
    orb_iter_input.WritePDB( out );
    fclose(out);
    
    assert(orb_iter_output.hasSelfEnergy());
    assert(orb_iter_output.hasQMEnergy());
    double energy___ex = 0.0;
    // EXTRACT & SAVE QM ENERGIES
    double energy___sf = orb_iter_output.getSelfEnergy();
    double energy_qmsf = orb_iter_output.getQMEnergy();
    double energy_qm__ = energy_qmsf - energy___sf ;
    thisIter->setQMSF(energy_qm__, energy___sf, energy___ex);
    _job->setEnergy_QMMM(thisIter->getQMEnergy(),thisIter->getGWBSEEnergy(), thisIter->getSFEnergy(),
                         thisIter->getQMMMEnergy());
    
    // EXTRACT & SAVE QMATOM DATA
    std::vector< ctp::QMAtom* > &atoms = *(orb_iter_output.getAtoms());
    
    thisIter->UpdatePosChrgFromQMAtoms(atoms, _job->getPolarTop()->QM0());

    LOG(ctp::logINFO,*_log) 
        << format("Summary - iteration %1$d:") % (iterCnt+1) << flush;
    LOG(ctp::logINFO,*_log)
        << format("... QM Size  = %1$d atoms") % int(atoms.size()) << flush;
    LOG(ctp::logINFO,*_log)
        << format("... E(QM)    = %1$+4.9e") % thisIter->getQMEnergy() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... E(GWBSE) = %1$+4.9e") % thisIter->getGWBSEEnergy() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... E(SF)    = %1$+4.9e") % thisIter->getSFEnergy() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... E(FM)    = %1$+4.9e") % thisIter->getFMEnergy() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... E(MM)    = %1$+4.9e") % thisIter->getMMEnergy() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... E(QMMM)  = %1$+4.9e") % thisIter->getQMMMEnergy() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... RMS(dR)  = %1$+4.9e") % thisIter->getRMSdR() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... RMS(dQ)  = %1$+4.9e") % thisIter->getRMSdQ() << flush;
    LOG(ctp::logINFO,*_log)
        << format("... SUM(dQ)  = %1$+4.9e") % thisIter->getSUMdQ() << flush;
    
    // CLEAN DIRECTORY
    _qmpack->CleanUp();

    
    /*
    int removed = boost::filesystem::remove_all(runFolder);
    if (removed > 0) 
        LOG(ctp::logDEBUG,*_log) << "Removed directory " << runFolder << flush;
    else 
        LOG(ctp::logWARNING,*_log) << "Could not remove dir " << runFolder << flush;
    */
    return 0;
     
}


template<class QMPackage>
QMMIter *QMAPEMachine<QMPackage>::CreateNewIter() {
    
    QMMIter *newIter = new QMMIter(_iters.size());
    this->_iters.push_back(newIter);
    return newIter;
}


template<class QMPackage>
bool QMAPEMachine<QMPackage>::EvaluateGWBSE(Orbitals &orb, string runFolder) {

	// for GW-BSE, we also need to parse the orbitals file
        
        _qmpack->ParseOrbitalsFile(&orb);
	//int _parse_orbitals_status = _qmpack->ParseOrbitalsFile(&orb);
	std::vector<int> _state_index;
        GWBSE _gwbse(&orb);
	_gwbse.Initialize( &_gwbse_options );
	if ( _state > 0 ){
	LOG(ctp::logDEBUG,*_log) << "Excited state via GWBSE: " <<  flush;
	LOG(ctp::logDEBUG,*_log) << "  --- type:              " << _type << flush;
	LOG(ctp::logDEBUG,*_log) << "  --- state:             " << _state << flush;
	if ( _has_osc_filter) { LOG(ctp::logDEBUG,*_log) << "  --- filter: osc.str. > " << _osc_threshold << flush; }
	if ( _has_dQ_filter) { LOG(ctp::logDEBUG,*_log) << "  --- filter: crg.trs. > " << _dQ_threshold << flush; }

	if ( _has_osc_filter && _has_dQ_filter ){
		LOG(ctp::logDEBUG,*_log) << "  --- WARNING: filtering for optically active CT transition - might not make sense... "  << flush;
	}

	// define own logger for GW-BSE that is written into a runFolder logfile
	ctp::Logger gwbse_logger(ctp::logDEBUG);
	gwbse_logger.setMultithreading(false);
	_gwbse.setLogger(&gwbse_logger);
	gwbse_logger.setPreface(ctp::logINFO,    (format("\nGWBSE INF ...") ).str());
	gwbse_logger.setPreface(ctp::logERROR,   (format("\nGWBSE ERR ...") ).str());
	gwbse_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...") ).str());
	gwbse_logger.setPreface(ctp::logDEBUG,   (format("\nGWBSE DBG ...") ).str());

	// actual GW-BSE run

	//bool _evaluate = _gwbse.Evaluate(&orb);
        _gwbse.Evaluate();

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
		const std::vector<ub::vector<double> >& TDipoles = orb.TransitionDipoles();
		for (unsigned _i=0; _i < TDipoles.size(); _i++ ) {

			double osc = (ub::inner_prod(TDipoles[_i],TDipoles[_i])) * 1.0 / 3.0 * (orb.BSESingletEnergies()(_i)) ;
			if ( osc > _osc_threshold ) _state_index.push_back(_i);
		}



	} else {

		if ( _type == "singlet" ){
		   for (unsigned _i=0; _i < orb.TransitionDipoles().size(); _i++ ) {
			   _state_index.push_back(_i);
		   }
		} else {
		   for (unsigned _i=0; _i < orb.BSETripletEnergies().size(); _i++ ) {
			   _state_index.push_back(_i);
		   }
		}
	}




        // filter according to charge transfer, go through list of excitations in _state_index
         if  (_has_dQ_filter ) {
            std::vector<int> _state_index_copy;
            if ( _type == "singlets" ){
             // go through list of singlets
            const std::vector<double>& dQ_fragA = orb.FragmentAChargesSingEXC();
            //const std::vector<double>& dQ_fragB = orb.FragmentBChargesSingEXC();
            for (unsigned _i=0; _i < _state_index.size(); _i++ ) {
                if ( std::abs(dQ_fragA[_i]) > _dQ_threshold ) {
                    _state_index_copy.push_back(_state_index[_i]);
                }
            } 
            _state_index = _state_index_copy;
            } else if ( _type == "triplets"){
              // go through list of triplets
            const std::vector<double>& dQ_fragA = orb.FragmentAChargesTripEXC();
            //const std::vector<double>& dQ_fragB = orb.FragmentBChargesTripEXC();
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
	

	} // only if state >0

	// calculate density matrix for this excited state
	ub::matrix<double> &_dft_orbitals = orb.MOCoefficients();
	// load DFT basis set (element-wise information) from xml file
	BasisSet dftbs;
	if ( orb.getDFTbasis() != "" ) {
	  dftbs.LoadBasisSet(orb.getDFTbasis());
	} else{
	dftbs.LoadBasisSet( _gwbse.get_dftbasis_name() );

	}
	LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " <<  orb.getDFTbasis()  << flush;





	// fill DFT AO basis by going through all atoms
	AOBasis dftbasis;
	dftbasis.AOBasisFill(&dftbs, orb.QMAtoms() );
	dftbasis.ReorderMOs(_dft_orbitals, orb.getQMpackage(), "xtp" );
	// TBD: Need to switch between singlets and triplets depending on _type
	ub::matrix<double> DMATGS=orb.DensityMatrixGroundState(_dft_orbitals);

	ub::matrix<double> DMAT_tot=DMATGS; // Ground state + hole_contribution + electron contribution

	if ( _state > 0 ){
	ub::matrix<real_gwbse>& BSECoefs = orb.BSESingletCoefficients();
	std::vector<ub::matrix<double> > DMAT = orb.DensityMatrixExcitedState( _dft_orbitals , BSECoefs, _state_index[_state-1]);
	DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1]; // Ground state + hole_contribution + electron contribution
	}

	// fill DFT AO basis by going through all atoms
	std::vector< ctp::QMAtom* >& Atomlist= orb.QMAtoms();

	Espfit esp=Espfit(_log);
        if (_run_gwbse){
        esp.setUseECPs(true);
        }
	esp.Fit2Density(Atomlist, DMAT_tot, dftbasis,dftbs,"medium");

	return true;
}


template<class QMPackage>
bool QMAPEMachine<QMPackage>::hasConverged() {
    
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
    
    LOG(ctp::logINFO,*_log) 
        << (format("Convergence check")) << flush;
    LOG(ctp::logINFO,*_log)
        << format("  o Converged dR ? %s") % (_convg_dR ? "True" : "False") << flush;
    LOG(ctp::logINFO,*_log) 
        << format("  o Converged dQ ? %s") % (_convg_dQ ? "True" : "False") << flush;
    LOG(ctp::logINFO,*_log) 
        << format("  o Converged QM ? %s") % (_convg_dE_QM ? "True" : "False") << flush;
    LOG(ctp::logINFO,*_log) 
        << format("  o Converged MM ? %s") % (_convg_dE_MM ? "True" : "False") << flush;
    
    return _isConverged;
}




// REGISTER QM PACKAGES
template class QMAPEMachine<QMPackage>;
    
    
    
}}
