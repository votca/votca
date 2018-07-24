/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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



#include <votca/xtp/qmapemachine.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/ctp/logger.h>

#include <votca/tools/elements.h>
#include <votca/xtp/espfit.h>

using boost::format;

namespace votca { namespace xtp {

QMAPEMachine::QMAPEMachine(ctp::XJob *job, ctp::Ewald3DnD *cape,
	 Property *opt, string sfx, int nst)
   : _job(job), _cape(cape),_isConverged(false) {
    
	// CONVERGENCE THRESHOLDS
    string key = sfx + ".qmmmconvg";
            _crit_dR =opt->ifExistsReturnElseReturnDefault<double>(key + ".dR",0.01);//nm
            _crit_dQ =opt->ifExistsReturnElseReturnDefault<double>(key + ".dQ",0.01);//e
            _crit_dE_QM =opt->ifExistsReturnElseReturnDefault<double>(key + ".dQdE_QM",0.001);//eV
            _crit_dE_MM  =opt->ifExistsReturnElseReturnDefault<double>(key + ".dE_MM",_crit_dE_QM);//eV
            _maxIter=opt->ifExistsReturnElseReturnDefault<int>(key + ".max_iter",32);

	// TASKS
	key = sfx + ".tasks";
		_run_ape = opt->ifExistsReturnElseReturnDefault<bool>(key+".run_ape",true);
		_run_dft = opt->ifExistsReturnElseReturnDefault<bool>(key+".run_dft",true);
		_run_gwbse = opt->ifExistsReturnElseReturnDefault<bool>(key+".run_gwbse",true);
        
        
        if(_run_dft){
            key=sfx+".dft";
            string dft_xml = opt->ifExistsReturnElseThrowRuntimeError<string>(key + ".dftengine");
            load_property_from_xml(_dft_options, dft_xml.c_str());
            _externalgridaccuracy=opt->ifExistsReturnElseReturnDefault<string>(key+".externalgrid","medium");
        }
                if(_run_gwbse){
	// GWBSE CONFIG
        key = sfx + ".gwbse";
		string gwbse_xml = opt->ifExistsReturnElseThrowRuntimeError<string>(key + ".gwbse_options");
		load_property_from_xml(_gwbse_options, gwbse_xml.c_str());

		_state = opt->ifExistsReturnElseReturnDefault<int>(key+".state",1);
                std::vector<string> choices={ "singlet","triplet"  };
		_type  = opt->ifExistsAndinListReturnElseThrowRuntimeError<string>(key+".type",choices);
		
                
		key = sfx + ".gwbse.filter";
                if (opt->exists(key + ".oscillator_strength") && _type != "triplet") {
                    _has_osc_filter = true;
                    _osc_threshold = opt->get(key + ".oscillator_strength").as<double> ();
                }
                if (opt->exists(key + ".localisation")) {
                    _has_loc_filter = true;
                    
                    string temp = opt->get(key + ".localisation").as<string> ();
                    Tokenizer tok_cleanup(temp, ", \n\t");
                    std::vector <std::string> strings_vec;
                    tok_cleanup.ToVector(strings_vec);
                    if (strings_vec.size()!=2){
                        throw runtime_error("qmmmachine: Fragment and localisation threshold are not separated");
                    }
                    if(strings_vec[0]=="a" || strings_vec[0]=="A"){
                        _localiseonA=true;
                    }else if(strings_vec[0]=="b" || strings_vec[0]=="B"){
                         _localiseonA=false;
                    }else{
                        throw runtime_error("qmmmachine: Fragment label not known, either A or B");
                    }
                    _loc_threshold=boost::lexical_cast<double>(strings_vec[1]);
                } 

                if (opt->exists(key + ".charge_transfer")) {
                    _has_dQ_filter = true;
                    _dQ_threshold = opt->get(key + ".charge_transfer").as<double> ();
                } 
                if(_has_dQ_filter && _has_loc_filter){
                    throw runtime_error("Cannot use localisation and charge_transfer filter at the same time.");
                }
                }
                return;
    }


QMAPEMachine::~QMAPEMachine() {
    
    std::vector<QMMIter*> ::iterator qit;
    for (qit = _iters.begin(); qit < _iters.end(); ++qit) {
        delete *qit;
    }
    _iters.clear();
    
    std::vector<ctp::PolarSeg*> ::iterator pit;
    for (pit = target_bg.begin(); pit < target_bg.end(); ++pit) {
        delete *pit;
    }
    target_bg.clear();
    
   
    for (pit = target_fg.begin(); pit < target_fg.end(); ++pit) {
        delete *pit;
    }
    target_fg.clear();
}


void QMAPEMachine::Evaluate(ctp::XJob *job) {
    
	// PREPARE JOB DIRECTORY
	string jobFolder = "qmapejob_" + boost::lexical_cast<string>(_job->getId())
					 + "_" + _job->getTag();
	bool created = boost::filesystem::create_directory(jobFolder);

	CTP_LOG(ctp::logINFO,*_log) << flush;
	if (created) {
		CTP_LOG(ctp::logINFO,*_log) << "Created directory " << jobFolder << flush;
        }

    CTP_LOG(ctp::logINFO,*_log)
       << format("... dR %1$1.4f dQ %2$1.4f QM %3$1.4f MM %4$1.4f IT %5$d")
       % _crit_dR % _crit_dQ % _crit_dE_QM % _crit_dE_MM % _maxIter << flush;
    
    // FIGURE OUT CHARGE + MULTIPLICITY
    double dQ = 0.0;
    for (unsigned i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
        dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
    }
    int chrg = round(dQ);
    int spin = ( (chrg < 0) ? -chrg:chrg ) % 2 + 1;
    CTP_LOG(ctp::logINFO,*_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;
    
    if(chrg!=0){
        throw runtime_error("Charged DFT calculations are not possible at the moment");
    }
    qminterface.setMultipoleSplitting(false,1);
    

    int iterCnt = 0;
    int iterMax = _maxIter;
    for ( ; iterCnt < iterMax; ++iterCnt) {
        CTP_LOG(ctp::logINFO,*_log) << "QMMM ITERATION:" << iterCnt+1<<" of "<<iterMax << flush;
        //bool code = Iterate(jobFolder, iterCnt);
        Iterate(jobFolder, iterCnt);
        if (hasConverged()) {
            CTP_LOG(ctp::logINFO,*_log) << "QMMM CONVERGED after:" << iterCnt+1<<" iterations."<< flush;
        
        break; }
    }
    
    if (iterCnt == iterMax-1 && !_isConverged) {
        CTP_LOG(ctp::logWARNING,*_log)
            << format("Not converged within %1$d iterations.") % iterMax;
    }
    
    return;
}


bool QMAPEMachine::Iterate(string jobFolder, int iterCnt) {

    // CREATE ITERATION OBJECT & SETUP RUN DIRECTORY
    QMMIter *thisIter = this->CreateNewIter();
    int iter = iterCnt;
    string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);
       
    CTP_LOG(ctp::logINFO,*_log) << flush;
    bool created = boost::filesystem::create_directory(runFolder);
    if (created) 
        CTP_LOG(ctp::logDEBUG,*_log) << "Created directory " << runFolder << flush;
    else
        CTP_LOG(ctp::logWARNING,*_log) << "Could not create directory " << runFolder << flush;
    
    Orbitals orb_iter_input;
    
    qminterface.GenerateQMAtomsFromPolarSegs(_job->getPolarTop(), orb_iter_input);
    
    
    
     
    if (iterCnt == 0) {
    dftengine.Initialize(_dft_options);
    dftengine.setLogger(_log);
    dftengine.ConfigureExternalGrid(_externalgridaccuracy);  
    dftengine.Prepare(orb_iter_input);
    SetupPolarSiteGrids(dftengine.getExternalGridpoints(),orb_iter_input.QMAtoms());
    }

    // COMPUTE POLARIZATION STATE WITH QM0(0)
    if (_run_ape) {
        if (iterCnt == 0) {
            _cape->ShowAgenda(_log);
            // Reset FGC, start from BGP state, apply FP fields (BG & FG)
            _cape->EvaluateInductionQMMM(true, true, true, true, true);
            // Add BG, do not add MM1 & QM0
            _cape->EvaluatePotential(target_bg, true, false, false);
        }
        // Do not add BG & QM0, add MM1
        _cape->EvaluatePotential(target_fg, false, true, false);
            }
    
    dftengine.setExternalGrid(ExtractElGrid_fromPolarsites(),ExtractNucGrid_fromPolarsites());
    
    if (_run_dft) {
    dftengine.Evaluate(orb_iter_input);
    }
  
	orb_iter_input.WriteXYZ(runFolder + "/Fullstructure.xyz","Full structure");

    // Run GWBSE
	if (_run_gwbse){
		EvaluateGWBSE(orb_iter_input, runFolder);
	}
        
       

	// COMPUTE POLARIZATION STATE WITH QM0(n+1)
	if (_run_ape) {
		// Update QM0 density: QM0(n) => QM0(n+1)
		// ...
        
        thisIter->UpdatePosChrgFromQMAtoms(orb_iter_input.QMAtoms(),_job->getPolarTop()->QM0());
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
}


QMMIter *QMAPEMachine::CreateNewIter() {
    
    QMMIter *newIter = new QMMIter(_iters.size());
    this->_iters.push_back(newIter);
    return newIter;
        }

bool QMAPEMachine::EvaluateGWBSE(Orbitals &orb, string runFolder) {

    std::vector<int> _state_index;
    GWBSE _gwbse(orb);
    _gwbse.setLogger(_log);
    
    if (_state > 0) {
        CTP_LOG(ctp::logDEBUG, *_log) << "Excited state via GWBSE: " << flush;
        CTP_LOG(ctp::logDEBUG, *_log) << "  --- type:              " << _type << flush;
        CTP_LOG(ctp::logDEBUG, *_log) << "  --- state:             " << _state << flush;
        if (_has_osc_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: osc.str. > " << _osc_threshold << flush;
        }
        if (_has_dQ_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: crg.trs. > " << _dQ_threshold << flush;
        }
        if (_has_loc_filter){
            if (_loc_threshold>0.5){
             CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: localisation > " << _loc_threshold << flush;
            }else{
                CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: localisation < " << _loc_threshold << flush;
            }
        }

        if (_has_osc_filter && _has_dQ_filter) {
            CTP_LOG(ctp::logDEBUG, *_log) << "  --- WARNING: filtering for optically active CT transition - might not make sense... " << flush;
        }

        // define own logger for GW-BSE that is written into a runFolder logfile
        ctp::Logger gwbse_logger(ctp::logDEBUG);
        gwbse_logger.setMultithreading(false);
        //gwbse_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
        gwbse_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
        gwbse_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
        gwbse_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());
        _gwbse.Initialize(_gwbse_options);
        // actual GW-BSE run

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


        if (_has_osc_filter) {

            // go through list of singlets
            const std::vector<double> oscs = orb.Oscillatorstrengths();
            for (unsigned _i = 0; _i < oscs.size(); _i++) {

                double osc =oscs[_i];
                if (osc > _osc_threshold) _state_index.push_back(_i);
            }
        } else {
            const VectorXfd& energies = (_type=="singlet") 
                        ? orb.BSESingletEnergies() : orb.BSETripletEnergies();
                       
                        for (unsigned _i = 0; _i < energies.size(); _i++) {
                            _state_index.push_back(_i);
                        }    
        }

       // filter according to charge transfer, go through list of excitations in _state_index
        if (_has_dQ_filter) {
            std::vector<int> _state_index_copy;
            const std::vector< Eigen::VectorXd >& dQ_frag= (_type=="singlet") 
            ? orb.getFragmentChargesSingEXC():orb.getFragmentChargesTripEXC();
            for (unsigned _i = 0; _i < _state_index.size(); _i++) {
                if (std::abs(dQ_frag[_state_index[_i]](0)) > _dQ_threshold) {
                    _state_index_copy.push_back(_state_index[_i]);
                }
            }
            _state_index = _state_index_copy;
        }
        else if (_has_loc_filter) {
            std::vector<int> _state_index_copy;
            const std::vector< Eigen::VectorXd >& popE= (_type=="singlet") 
            ? orb.getFragment_E_localisation_singlet():orb.getFragment_E_localisation_triplet();
            const std::vector< Eigen::VectorXd >& popH= (_type=="singlet") 
            ? orb.getFragment_H_localisation_singlet():orb.getFragment_H_localisation_triplet();
            if(_localiseonA){
                for (unsigned _i = 0; _i < _state_index.size(); _i++) {
                    if (popE[_state_index[_i]](0) > _loc_threshold && popH[_state_index[_i]](0) > _loc_threshold ) {
                        _state_index_copy.push_back(_state_index[_i]);
                    }
                }
            }else{
                for (unsigned _i = 0; _i < _state_index.size(); _i++) {
                    if (popE[_state_index[_i]](1) > _loc_threshold && popH[_state_index[_i]](1) > _loc_threshold ) {
                        _state_index_copy.push_back(_state_index[_i]);
                    }
                }
            }
            _state_index = _state_index_copy;
        }
        if (_state_index.size() < 1) {
            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " WARNING: FILTER yielded no state. Taking lowest excitation"<< flush;
            _state_index.push_back(0);
        }else{
                        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded state"<<_type<<":"<<_state_index[_state - 1]+1<< flush;
        }
    } // only if state >0
    
   
    
    // load DFT basis set (element-wise information) from xml file
    BasisSet dftbs;
    
    dftbs.LoadBasisSet(orb.getDFTbasis());
    
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb.getDFTbasis() << flush;

    // fill DFT AO basis by going through all atoms
    AOBasis dftbasis;
    dftbasis.AOBasisFill(dftbs, orb.QMAtoms() );
    // TBD: Need to switch between singlets and triplets depending on _type
   Eigen::MatrixXd DMATGS = orb.DensityMatrixGroundState();
    Eigen::MatrixXd DMAT_tot = DMATGS; // Ground state + hole_contribution + electron contribution

    if (_state > 0) {
        
        std::vector< Eigen::MatrixXd > DMAT = orb.DensityMatrixExcitedState(_type, _state_index[_state - 1]);
        DMAT_tot = DMAT_tot - DMAT[0] + DMAT[1]; // Ground state + hole_contribution + electron contribution
    }

    // fill DFT AO basis by going through all atoms
    
    Espfit esp = Espfit(_log);
    esp.Fit2Density(orb.QMAtoms(), DMAT_tot, dftbasis,"medium");
    return true;
}


void QMAPEMachine::SetupPolarSiteGrids(const std::vector<const vec *>& gridpoints,const std::vector<QMAtom*>& atoms){
    NumberofAtoms=0;
    std::vector<QMAtom*>::const_iterator qmt;
    std::vector<ctp::APolarSite*> sites1;
    std::vector<ctp::APolarSite*> sites2;
    
    for(qmt=atoms.begin();qmt!=atoms.end();++qmt){
        NumberofAtoms++;
        sites1.push_back(qminterface.Convert(*qmt));
        sites2.push_back(qminterface.Convert(*qmt));
    }
    
    
    std::vector<const vec *>::const_iterator grt;
    for (grt=gridpoints.begin();grt<gridpoints.end();++grt){
        ctp::APolarSite* site1=new ctp::APolarSite();
        ctp::APolarSite* site2=new ctp::APolarSite();
        tools::vec pos=*(*grt)*tools::conv::bohr2nm;
        site1->setPos(pos);
        site2->setPos(pos);
        sites1.push_back(site1);
        sites2.push_back(site2);
    }
    
    
    target_bg.push_back( new ctp::PolarSeg(0, sites1));
    target_fg.push_back( new ctp::PolarSeg(0, sites2));
    
    return;
        }

std::vector<double> QMAPEMachine::ExtractNucGrid_fromPolarsites() {
    std::vector<double> gridpoints;
    double int2hrt = tools::conv::int2V * tools::conv::ev2hrt;
    for (unsigned i = 0; i < target_bg.size(); ++i) {
        for (unsigned j = 0; j < NumberofAtoms; ++j) {
            ctp::PolarSeg* seg1 = target_bg[i];
            ctp::PolarSeg* seg2 = target_fg[i];
            ctp::APolarSite* site1 = (*seg1)[j];
            ctp::APolarSite* site2 = (*seg2)[j];
            double value = (site1->getPhi() + site2->getPhi()) * int2hrt;
            gridpoints.push_back(value);
        }
    }
    return gridpoints;
}

std::vector<double> QMAPEMachine::ExtractElGrid_fromPolarsites() {
    std::vector<double> gridpoints;
    double int2hrt = tools::conv::int2V * tools::conv::ev2hrt;
    for (unsigned i = 0; i < target_bg.size(); ++i) {
        for (unsigned j = NumberofAtoms; j < target_bg[i]->size(); ++j) {
            ctp::PolarSeg* seg1 = target_bg[i];
            ctp::PolarSeg* seg2 = target_fg[i];
            ctp::APolarSite* site1 = (*seg1)[j];
            ctp::APolarSite* site2 = (*seg2)[j];
            double value = (site1->getPhi() + site2->getPhi()) * int2hrt;
            gridpoints.push_back(value);
        }
    }
    return gridpoints;
}




bool QMAPEMachine::hasConverged() {
    
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
    
    CTP_LOG(ctp::logINFO,*_log) 
        << (format("Convergence check")) << flush;
    CTP_LOG(ctp::logINFO,*_log)
        << format("  o Converged dR ? %s") % (_convg_dR ? "True" : "False") << flush;
    CTP_LOG(ctp::logINFO,*_log) 
        << format("  o Converged dQ ? %s") % (_convg_dQ ? "True" : "False") << flush;
    CTP_LOG(ctp::logINFO,*_log) 
        << format("  o Converged QM ? %s") % (_convg_dE_QM ? "True" : "False") << flush;
    CTP_LOG(ctp::logINFO,*_log) 
        << format("  o Converged MM ? %s") % (_convg_dE_MM ? "True" : "False") << flush;
    
    return _isConverged;
}




    
    
    
}}
