// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/qmapemachine.h>
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
QMAPEMachine<QMPackage>::QMAPEMachine(XJob *job, Ewald3DnD *cape, QMPackage *qmpack,
	 Property *opt, string sfx, int nst)
   : _subthreads(nst),_job(job), _qmpack(qmpack), _cape(cape), 
	  _grid_fg(true,true,true), _grid_bg(true,true,true),
     _fitted_charges(true,true,true),_isConverged(false) {
    
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
    
    vector<QMAPEIter*> ::iterator qit;
    for (qit = _iters.begin(); qit < _iters.end(); ++qit) {
        delete *qit;
    }
    _iters.clear();
}


template<class QMPackage>
void QMAPEMachine<QMPackage>::Evaluate(XJob *job) {
    
	// PREPARE JOB DIRECTORY
	string jobFolder = "job_" + boost::lexical_cast<string>(_job->getId())
					 + "_" + _job->getTag();
	bool created = boost::filesystem::create_directory(jobFolder);

	LOG(logINFO,*_log) << flush;
	if (created) {
		LOG(logINFO,*_log) << "Created directory " << jobFolder << flush;
        }

    LOG(logINFO,*_log)
       << format("... dR %1$1.4f dQ %2$1.4f QM %3$1.4f MM %4$1.4f IT %5$d")
       % _crit_dR % _crit_dQ % _crit_dE_QM % _crit_dE_MM % _maxIter << flush;
    
    // FIGURE OUT CHARGE + MULTIPLICITY
    double dQ = 0.0;
    for (unsigned i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
        dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
    }
    int chrg = round(dQ);
    int spin = ( (chrg < 0) ? -chrg:chrg ) % 2 + 1;
    LOG(logINFO,*_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;

    // SET ITERATION-TIME CONSTANTS
    _qmpack->setCharge(chrg);
    _qmpack->setSpin(spin);

    // GENERATE GRIDS
    // ... TODO ...
    // Generate QM atoms from _job->getPolarTop()->QM0();
    // Move Iter::GenerateQMAtomsFromPolarSegs to QMMachine
    // Generate grids, store as member
    Orbitals basisforgrid;
    vector<PolarSeg*> dummy;
    
    GenerateQMAtomsFromPolarSegs(_job->getPolarTop()->QM0(),dummy,basisforgrid);
    
    
    _grid_bg = Grid(true,true,true);   
    _grid_bg.setAtomlist(&basisforgrid.QMAtoms());
    _grid_bg.setCutoffshifts(0,-0.5);    
    _grid_bg.setSpacing(0.2);   
    _grid_bg.setCubegrid(true);   
    _grid_bg.setupgrid();
    
    LOG(logINFO,*_log) << "Created internal background grid with " << _grid_bg.getsize() <<" points."<< flush;

    _grid_fg=_grid_bg;
    
    LOG(logINFO,*_log) << "Created internal foreground grid with " << _grid_fg.getsize() <<" points."<< flush;
    

            
    _fitted_charges = Grid(true,false,false);
    _fitted_charges.setAtomlist(&basisforgrid.QMAtoms());
    _fitted_charges.setCutoffs(7,0);
    _fitted_charges.setupradialgrid(1);
    //_fitted_charges.setCutoffshifts(8,2);
    //_fitted_charges.setSpacing(3);
    
    //_fitted_charges.setupgrid();
    
    LOG(logINFO,*_log) << "Created " << _fitted_charges.getsize() <<" charge positions."<< flush;
    LOG(logINFO,*_log) << flush;
    _exportgridtofile=true;
    if (_exportgridtofile){
    _grid_bg.printGridtoxyzfile("grid.xyz");
    _fitted_charges.printGridtoxyzfile("grid2.xyz");
    }
    
    int iterCnt = 0;
    int iterMax = _maxIter;
    for ( ; iterCnt < iterMax; ++iterCnt) {
        
        //bool code = Iterate(jobFolder, iterCnt);
        Iterate(jobFolder, iterCnt);
        if (hasConverged()) { break; }
    }
    
    if (iterCnt == iterMax-1 && !_isConverged) {
        LOG(logWARNING,*_log)
            << format("Not converged within %1$d iterations.") % iterMax;
    }
    
    return;
}


template<class QMPackage>
bool QMAPEMachine<QMPackage>::Iterate(string jobFolder, int iterCnt) {

    // CREATE ITERATION OBJECT & SETUP RUN DIRECTORY
    QMAPEIter *thisIter = this->CreateNewIter();
    int iter = iterCnt;
    string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);
       
    LOG(logINFO,*_log) << flush;
    bool created = boost::filesystem::create_directory(runFolder);
    if (created) 
        LOG(logDEBUG,*_log) << "Created directory " << runFolder << flush;
    else
        LOG(logWARNING,*_log) << "Could not create directory " << runFolder << flush;

    // COMPUTE POLARIZATION STATE WITH QM0(0)
    if (_run_ape) {
		if (iterCnt == 0) {
			_cape->ShowAgenda(_log);
			// Reset FGC, start from BGP state, apply FP fields (BG & FG)
			//_cape->EvaluateInductionQMMM(true, true, true, true, true);
		}
    
     
        //vec pos1=0.5*(vec(_fitted_charges.getGrid()[0])+vec(_fitted_charges.getGrid()[1]));
        //vec pos2=0.5*(vec(_fitted_charges.getGrid()[10])+vec(_fitted_charges.getGrid()[11]));
        double q1=-1.0;
        //double q2=1.0;
     vec pos=vec(45.135,2.484,-5.54117);
		
       
        vector<APolarSite*>::iterator pit;
        for (pit=_grid_bg.Sites().begin();pit!=_grid_bg.Sites().end();++pit){
            //double dist1=abs((*pit)->getPos()-pos1);
            //double dist2=abs((*pit)->getPos()-pos2);
            double dist=abs((*pit)->getPos()-pos);
            //double potential=q1/dist1+q2/dist2;  
            double potential=q1/dist;
            (*pit)->setPhi(potential,0.0);
        }         
        
        
		vector< PolarSeg* > target_bg;     
        target_bg.push_back(_grid_bg.getSeg());
        
        vector< PolarSeg* > target_fg;
        target_fg.push_back(_grid_fg.getSeg());
        
       
		if (iterCnt == 0) {
			// Add BG, do not add MM1 & QM0
			//_cape->EvaluatePotential(target_bg, true, false, false);
		}
		// Do not add BG & QM0, add MM1
		//_cape->EvaluatePotential(target_fg, false, true, false);
    }
    
    
      /*
        _grid_bg.getSeg()->WriteMPS("test_bg.mps", "TEST");
        cout << endl << "Done bg. " << endl;
        
        cout << endl << "Done ... " << endl;
        _grid_fg.getSeg()->WriteMPS("test_fg.mps", "TEST");
        cout << endl << "Done fg. " << endl;
        
         cout << endl << "Done ... " << endl;
        _fitted_charges.getSeg()->WriteMPS("test_charges.mps", "TEST");
        cout << endl << "Done charges. " << endl;
    */  
        
    //_grid_fg.readgridfromCubeFile("cubefile_fg.cub");
    //_grid_bg.readgridfromCubeFile("cubefile_bg.cub");
    _grid_fg.printgridtoCubefile("cubefile_fg.cub");
    _grid_bg.printgridtoCubefile("cubefile_bg.cub");
    Grid test;
    test.readgridfromCubeFile("cubefile_fg.cub");
    test.printgridtoCubefile("cubefile_fg_out.cub");

    // COMPUTE WAVEFUNCTION & QM ENERGY
    // Generate charge shell from potentials
    vector<PolarSeg*> &qm =_job->getPolarTop()->QM0();
    vector<PolarSeg*> mm_fitted;
    Espfit fitcharges;
    fitcharges.setLog(_log);
    double netchargefit=0.0;

    fitcharges.FitAPECharges(_grid_bg,_grid_fg,_fitted_charges,netchargefit);
    mm_fitted.push_back(_fitted_charges.getSeg());

    
    vector<PolarSeg*> dummy;
    Orbitals basisforgrid;
    GenerateQMAtomsFromPolarSegs(_job->getPolarTop()->QM0(),dummy,basisforgrid);
    Grid visgrid_fit=Grid(_grid_bg);
  
    fitcharges.EvaluateAPECharges(visgrid_fit,_fitted_charges);

    visgrid_fit.printgridtoCubefile("cubefile_fit.cub");

  
    
  
    exit(0);
       
    
    
    // Run DFT
    Orbitals orb_iter_input;
    vector<Segment*> empty;
    GenerateQMAtomsFromPolarSegs(qm, mm_fitted, orb_iter_input);
   
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
    //_cape->ShowAgenda(_log);
    //_cape->EvaluateInductionQMMM(true, true, true, true, true);
    //_cape->EvaluateEnergyQMMM();
    //_cape->ShowEnergySplitting(_log);

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
    
    LOG(logDEBUG,*_log) << "Writing input file " << runFolder << flush;
    
    _qmpack->WriteInputFile(empty, &orb_iter_input);
         
    // RUN HERE (OVERRIDE - COPY EXISTING LOG-FILE)
    //string cpstr = "cp e_1_n.log " + path_logFile;
    //int sig = system(cpstr.c_str());
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
QMAPEIter *QMAPEMachine<QMPackage>::CreateNewIter() {
    
    QMAPEIter *newIter = new QMAPEIter(_iters.size());
    this->_iters.push_back(newIter);
    return newIter;
}


template<class QMPackage>
bool QMAPEMachine<QMPackage>::EvaluateGWBSE(Orbitals &orb, string runFolder) {

	// for GW-BSE, we also need to parse the orbitals file
        _qmpack->ParseOrbitalsFile(&orb);
	//int _parse_orbitals_status = _qmpack->ParseOrbitalsFile(&orb);
	std::vector<int> _state_index;
	_gwbse.Initialize( &_gwbse_options );
	if ( _state > 0 ){
	LOG(logDEBUG,*_log) << "Excited state via GWBSE: " <<  flush;
	LOG(logDEBUG,*_log) << "  --- type:              " << _type << flush;
	LOG(logDEBUG,*_log) << "  --- state:             " << _state << flush;
	if ( _has_osc_filter) { LOG(logDEBUG,*_log) << "  --- filter: osc.str. > " << _osc_threshold << flush; }
	if ( _has_dQ_filter) { LOG(logDEBUG,*_log) << "  --- filter: crg.trs. > " << _dQ_threshold << flush; }

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

	//bool _evaluate = _gwbse.Evaluate(&orb);
        _gwbse.Evaluate(&orb);

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
		const std::vector<std::vector<double> >& TDipoles = orb.TransitionDipoles();
		for (unsigned _i=0; _i < TDipoles.size(); _i++ ) {

			double osc = (TDipoles[_i][0] * TDipoles[_i][0] + TDipoles[_i][1] * TDipoles[_i][1] + TDipoles[_i][2] * TDipoles[_i][2]) * 1.0 / 3.0 * (orb.BSESingletEnergies()[_i]) ;
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

 if ( _has_osc_filter ){
            
            // go through list of singlets
            const std::vector<std::vector<double> >& TDipoles = orb.TransitionDipoles();
            for (unsigned _i=0; _i < TDipoles.size(); _i++ ) {
                
                double osc = (TDipoles[_i][0] * TDipoles[_i][0] + TDipoles[_i][1] * TDipoles[_i][1] + TDipoles[_i][2] * TDipoles[_i][2]) * 1.0 / 3.0 * (orb.BSESingletEnergies()[_i]) ;
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
	// - output its energy
        /*
	double energy___ex = 0.0;
	if ( _type == "singlet" ){
		energy___ex = orb.BSESingletEnergies()[_state_index[_state-1]]*13.6058; // to eV
	} else if ( _type == "triplet" ) {
		energy___ex = orb.BSETripletEnergies()[_state_index[_state-1]]*13.6058; // to eV
	}
*/
	// ub::matrix<double> &_dft_orbitals_GS = orb_iter_output.MOCoefficients();
	// int _parse_orbitals_status_GS = _qmpack->ParseOrbitalsFile( &orb_iter_output );

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
	LOG(logDEBUG, *_log) << TimeStamp() << " Loaded DFT Basis Set " <<  orb.getDFTbasis()  << flush;





	// fill DFT AO basis by going through all atoms
	AOBasis dftbasis;
	dftbasis.AOBasisFill(&dftbs, orb.QMAtoms() );
	dftbasis.ReorderMOs(_dft_orbitals, orb.getQMpackage(), "votca" );
	// TBD: Need to switch between singlets and triplets depending on _type
	ub::matrix<double> &DMATGS=orb.DensityMatrixGroundState(_dft_orbitals);

	ub::matrix<double> DMAT_tot=DMATGS; // Ground state + hole_contribution + electron contribution

	if ( _state > 0 ){
	ub::matrix<float>& BSECoefs = orb.BSESingletCoefficients();
	std::vector<ub::matrix<double> > &DMAT = orb.DensityMatrixExcitedState( _dft_orbitals , BSECoefs, _state_index[_state-1]);
	DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1]; // Ground state + hole_contribution + electron contribution
	}

	// fill DFT AO basis by going through all atoms
	vector< QMAtom* >& Atomlist= orb.QMAtoms();

	Espfit esp;

	esp.setLog(_log);
	esp.Fit2Density(Atomlist, DMAT_tot, dftbasis,dftbs,"medium",false);

	return true;
}


template<class QMPackage>
bool QMAPEMachine<QMPackage>::hasConverged() {
    
    _convg_dR = false;
    _convg_dQ = false;
    _convg_dE_QM = false;
    _convg_dE_MM = false;
    
    if (_iters.size() > 1) {
        
        QMAPEIter *iter_0 = _iters[_iters.size()-2];
        QMAPEIter *iter_1 = _iters[_iters.size()-1];
        
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
        << (format("Convergence check")) << flush;
    LOG(logINFO,*_log)
        << format("  o Converged dR ? %s") % (_convg_dR ? "True" : "False") << flush;
    LOG(logINFO,*_log) 
        << format("  o Converged dQ ? %s") % (_convg_dQ ? "True" : "False") << flush;
    LOG(logINFO,*_log) 
        << format("  o Converged QM ? %s") % (_convg_dE_QM ? "True" : "False") << flush;
    LOG(logINFO,*_log) 
        << format("  o Converged MM ? %s") % (_convg_dE_MM ? "True" : "False") << flush;
    
    return _isConverged;
}


void QMAPEIter::ConvertPSitesToQMAtoms(vector< PolarSeg* > &psegs,
                                       vector< QMAtom * > &qmatoms) {
    
    assert(qmatoms.size() == 0);    
    return;   
}


void QMAPEIter::ConvertQMAtomsToPSites(vector< QMAtom* > &qmatoms,
                                       vector< PolarSeg* > &psegs) {
    assert(qmatoms.size() == 0);
    return;
}


void QMAPEIter::UpdatePosChrgFromQMAtoms(vector< QMAtom* > &qmatoms,
                                         vector< PolarSeg* > &psegs) {
    
    double AA_to_NM = 0.1; // Angstrom to nanometer
    
    double dR_RMS = 0.0;
    double dQ_RMS = 0.0;
    double dQ_SUM = 0.0;
    
    for (unsigned i = 0, qac = 0; i < psegs.size(); ++i) {
        PolarSeg *pseg = psegs[i];
        for (unsigned j = 0; j < pseg->size(); ++j, ++qac) {
            
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

template<class QMPackage>
void QMAPEMachine<QMPackage>::GenerateQMAtomsFromPolarSegs(vector<PolarSeg*> &qm,
	vector<PolarSeg*> &mm, Orbitals &orb) {
    
    double AA_to_NM = 0.1; // Angstrom to nanometer
    
    // QM REGION
    for (unsigned i = 0; i < qm.size(); ++i) {
        vector<APolarSite*> *pseg = qm[i];
        for (unsigned j = 0; j < pseg->size(); ++j) {
            APolarSite *aps = (*pseg)[j];
            string type = "qm";
            vec pos = aps->getPos()/AA_to_NM;
            double Q = 0.0;
            orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, false);
        }
    }
    
    // MM REGION (EXPANDED VIA PARTIAL CHARGES)
    for (unsigned i = 0; i < mm.size(); ++i) {
    	vector<APolarSite*> *pseg = mm[i];
        for (unsigned j = 0; j < pseg->size(); ++j) {
            APolarSite *aps = (*pseg)[j];
            string type = "mm";
            vec pos = aps->getPos()/AA_to_NM;
            double Q = aps->getQ00();
            orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, true);
        }
    }
    
    return;
}\





void QMAPEIter::setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM) {
    
    _hasdRdQ = true;    
    _dR_RMS = dR_RMS;
    _dQ_RMS = dQ_RMS;
    _dQ_SUM = dQ_SUM;
    return;
}


void QMAPEIter::setQMSF(double energy_QM, double energy_SF, double energy_GWBSE) {
    
    _hasQM = true;
    _e_QM = energy_QM;
    _e_SF = energy_SF;    

    _hasGWBSE = true;
    _e_GWBSE = energy_GWBSE;
   
    return;
}


void QMAPEIter::setE_FM(double ef00, double ef01, double ef02, 
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


double QMAPEIter::getMMEnergy() {
    
    assert(_hasMM);
    return _ef_11 + _ef_12 + _em_1_ + _em_2_;
}


double QMAPEIter::getQMMMEnergy() {
    
    assert(_hasQM && _hasMM && _hasGWBSE);    
    return _e_QM + + _e_GWBSE + _ef_11 + _ef_12 + _em_1_ + _em_2_;    
}


// REGISTER QM PACKAGES
template class QMAPEMachine<QMPackage>;
    
    
    
}}
