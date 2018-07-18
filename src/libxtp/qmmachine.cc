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


#include <votca/xtp/qmmachine.h>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/elements.h>
#include <votca/xtp/espfit.h>


using boost::format;

namespace votca {
    namespace xtp {

        template<class QMPackage>
        QMMachine<QMPackage>::QMMachine(ctp::XJob *job, ctp::XInductor *xind, QMPackage *qmpack,
                Property *opt, string sfx, int nst, bool mav)
        : _job(job), _xind(xind), _qmpack(qmpack), _subthreads(nst),
        _isConverged(false) {

            string key = sfx + ".qmmmconvg";

            _crit_dR = opt->ifExistsReturnElseReturnDefault<double>(key + ".dR", 0.01); //nm
            _crit_dQ = opt->ifExistsReturnElseReturnDefault<double>(key + ".dQ", 0.01); //e
            _crit_dE_QM = opt->ifExistsReturnElseReturnDefault<double>(key + ".dQdE_QM", 0.001); //eV
            _crit_dE_MM = opt->ifExistsReturnElseReturnDefault<double>(key + ".dE_MM", _crit_dE_QM); //eV
            _maxIter = opt->ifExistsReturnElseReturnDefault<int>(key + ".max_iter", 32);
            
             _alpha = opt->ifExistsReturnElseReturnDefault<double>(key + ".mixing", 0.0); 
            
           
  


            key = sfx;
            bool split_dpl = opt->ifExistsReturnElseReturnDefault<bool>(key + ".split_dpl", true);
            double dpl_spacing = opt->ifExistsReturnElseReturnDefault<double>(key + ".dpl_spacing", 1e-3);
            qminterface.setMultipoleSplitting(split_dpl, dpl_spacing);

            // check for archiving
            std::string _archiving_string = opt->ifExistsReturnElseReturnDefault<std::string>(key + ".archiving", "");
            if (_archiving_string.find("iterations") != std::string::npos) {
                _do_archive = true;
            } else {
                _do_archive = false;
            }

            // check for static or polarized qmmm
            key = sfx + ".tholemodel";
            _static_qmmm = true;
             qmpack->setWithPolarization(false);
            if (opt->exists(key + ".induce")) {
                bool induce = opt->get(key + ".induce").as<bool>();
                
                if(induce){
                  qmpack->setWithPolarization(true);
                  qmpack->setDipoleSpacing(dpl_spacing);
                }
           
                cout<<"STATIC "<<induce<<endl;
                _static_qmmm = !induce;
            }
           


            // GDMA options
            key = sfx + ".gdma";
            if (opt->exists(key)) {
                _do_gdma = true;
                string _gdma_xml = opt->get(key).as<string> ();
                //cout << endl << "... ... Parsing " << _package_xml << endl ;
                load_property_from_xml(_gdma_options, _gdma_xml.c_str());
            } else {
                _do_gdma = false;
            }


            // GWBSE options
            key = sfx + ".gwbse";
            if (opt->exists(key)) {
                _do_gwbse = true;

            // PUT IN some useful type definitions for excited state
            // - singlet -> singlet exciton
            // - triplet -> triplet exciton
            // - quasiparticle -> GW quasiparticle (diagqp)
            // - kohn-sham -> DFT MO

                string _gwbse_xml = opt->get(key + ".gwbse_options").as<string> ();
                //cout << endl << "... ... Parsing " << _package_xml << endl ;
                load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());
                _state = opt->get(key + ".state").as< int >();
                _type = opt->get(key + ".type").as< string >();
                if (_type != "singlet" && _type != "triplet" && _type != "quasiparticle") {
                    throw runtime_error(" Invalid excited state type! " + _type);
                }

                key = sfx + ".gwbse.filter";
                if (opt->exists(key + ".oscillator_strength") && _type != "triplet") {
                    _has_osc_filter = true;
                    _osc_threshold = opt->get(key + ".oscillator_strength").as<double> ();
                }
                if (opt->exists(key + ".overlap")) {
                    _has_overlap_filter = true;
                    // _osc_threshold = opt->get(key + ".oscillator_strength").as<double> ();
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
            } else {
                _do_gwbse = false;
            }

            return;
        }

        template<class QMPackage>
        QMMachine<QMPackage>::~QMMachine() {

            std::vector<QMMIter*> ::iterator qit;
            for (qit = _iters.begin(); qit < _iters.end(); ++qit) {
                delete *qit;
            }
            _iters.clear();
        }

        template<class QMPackage>
        int QMMachine<QMPackage>::Evaluate(ctp::XJob *job) {

            CTP_LOG(ctp::logINFO, *_log)
                    << format("... dR %1$1.4f dQ %2$1.4f QM %3$1.4f MM %4$1.4f IT %5$d")
                    % _crit_dR % _crit_dQ % _crit_dE_QM % _crit_dE_MM % _maxIter << flush;

            // FIGURE OUT CHARGE + MULTIPLICITY
            double dQ = 0.0;
            for (unsigned int i = 0; i < _job->getPolarTop()->QM0().size(); ++i) {
                dQ += _job->getPolarTop()->QM0()[i]->CalcTotQ();
            }
            int chrg = round(dQ);
            int spin = ((chrg < 0) ? -chrg : chrg) % 2 + 1;
            CTP_LOG(ctp::logINFO, *_log) << "... Q = " << chrg << ", 2S+1 = " << spin << flush;


            // PREPARE JOB DIRECTORY
            string jobFolder = "xjob_" + boost::lexical_cast<string>(_job->getId())
                    + "_" + _job->getTag();
            bool created = boost::filesystem::create_directory(jobFolder);
            if (created) {
                CTP_LOG(ctp::logINFO, *_log) << "Created directory " << jobFolder << flush;
            }


            // SET ITERATION-TIME CONSTANTS
            // TO ADJUST

            _qmpack->setCharge(chrg);
            _qmpack->setSpin(spin);


            int iterCnt = 0;
            int iterMax = _maxIter;
            for (; iterCnt < iterMax; ++iterCnt) {

                // check for polarized QM/MM convergence
            int info= Iterate(jobFolder, iterCnt);
            if(info!=0){
                 CTP_LOG(ctp::logERROR, *_log)
                        << format("Iterating job failed!")<<flush;
                 return 1;
            }
                if (_static_qmmm) {
                    _isConverged = true;
                    break;
                } else if (hasConverged()) {
                    break;
                }

            }

            if (iterCnt == iterMax - 1 && !_isConverged) {
                CTP_LOG(ctp::logWARNING, *_log)
                        << format("Not converged within %1$d iterations.") % iterMax;
                return 2;
            }

            return 0;
        }

        template<class QMPackage>
        bool QMMachine<QMPackage>::Iterate(string jobFolder, int iterCnt) {

            // CREATE ITERATION OBJECT & SETUP RUN DIRECTORY
            QMMIter *thisIter = this->CreateNewIter();
            int iter = iterCnt;
            string runFolder = jobFolder + "/iter_" + boost::lexical_cast<string>(iter);

            bool created = boost::filesystem::create_directory(runFolder);
            if (created)
                CTP_LOG(ctp::logDEBUG, *_log) << "Created directory " << runFolder << flush;
            else
                CTP_LOG(ctp::logWARNING, *_log) << "Could not create directory " << runFolder << flush;


            // RUN CLASSICAL INDUCTION & SAVE
            _job->getPolarTop()->PrintPDB(runFolder + "/QM0_MM1_MM2.pdb");
            _xind->Evaluate(_job);



            assert(_xind->hasConverged());
            thisIter->setE_FM(_job->getEF00(), _job->getEF01(), _job->getEF02(),
                    _job->getEF11(), _job->getEF12(), _job->getEM0(),
                    _job->getEM1(), _job->getEM2(), _job->getETOT());


            std::vector<ctp::Segment*> empty;
            /* Translate atoms in QM0() to QMAtoms in orbitals object
             * to be used in writing the QM input files for the
             * external QMPackages or directly in internal DFT engine.
             * DRAWBACK: QMAtom positions are fixed for all iterations
             *           unless the UPDATE function at the end of the
             *           iteration updates orb_iter_input (TO BE TESTED)
             */
            if (iterCnt == 0) qminterface.GenerateQMAtomsFromPolarSegs(_job->getPolarTop(), orb_iter_input);

            /* Generate list of polar segments in the MM1() and MM2()
             * region to be used in writing the background multipoles
             * in the external QMPackages or in the direct evaluation
             * in the internal DFT engine.
             */
            std::vector<ctp::PolarSeg*> MultipolesBackground = qminterface.GenerateMultipoleList( _job->getPolarTop() );

            _qmpack->setMultipoleBackground( MultipolesBackground );

            // setting RUNDIR for the external QMPackages, dummy for internal
            _qmpack->setRunDir(runFolder);

            /* Call to WriteInputFile writes the appropriate input files
             * for the respective external QMPackages. For the internal
             * DFT engine, this function sets logger and runs DFTENGINE's
             * Prepare() function. ONLY the first iteration, this will
             * initialize the atoms, basissets, ecps, etc. In all
             * subsequent iterations it will recalculate all the "static"
             * AOmatrices (overlap, kinetic energy, nuc/ecp) which is
             * strictly unnecessary but at the same time also those
             * for external point charges, dipoles, and quadrupoles.
             * Since in polarized calculations, the dipoles can change
             * this recalculation is required. Should be split off the
             * Prepare() function for efficiency.
             */
            CTP_LOG(ctp::logDEBUG, *_log) << "Writing input file " << runFolder << flush;
            _qmpack->WriteInputFile(empty, &orb_iter_input);


            FILE *out;
            out = fopen((runFolder + "/system.pdb").c_str(), "w");
            orb_iter_input.WritePDB(out);
            fclose(out);

            /* Runs the external QMPackage or the self-consistent part of
             * DFTENGINE
             */
            _qmpack->Run( &orb_iter_input );

            /* Parse information from the LOGFILE into orbitals, if
             * external QMPackage is run. Dummy for internal DFTENGINE.
             * The QMPackage's MOcoefficients are not automatically
             * parsed in DFT-only calculations and ESP fits for
             * polarized QMMM are expected in the LOGFILE.
             * IF the internal ESPFITs should be used, the MOcoefficients
             * need to be parsed too.
             */
            bool success=_qmpack->ParseLogFile(orb_iter_input);
            if(!success){
                return 1;
            }
            
            // GW-BSE starts here
            double energy___ex = 0.0;
            std::vector<int> _state_index;

            if (_do_gwbse) {

                /* Parses the MOcoefficients from the external QMPackages
                 * for GW-BSE. Internal DFTENGINE has stored coefficients
                 * into orb_iter_input already, so this is a dummy for that.
                 */
                _qmpack->ParseOrbitalsFile(orb_iter_input);
                orb_iter_input.setDFTbasis(_qmpack->getBasisSetName());

                // Get a GWBSE object
                GWBSE _gwbse = GWBSE(orb_iter_input);
                // define own logger for GW-BSE that is written into a runFolder logfile
                ctp::Logger gwbse_logger(ctp::logDEBUG);
                gwbse_logger.setMultithreading(false);
                //_gwbse.setLogger(_log);
                _gwbse.setLogger(&gwbse_logger);
                gwbse_logger.setPreface(ctp::logINFO, (format("\nGWBSE INF ...")).str());
                gwbse_logger.setPreface(ctp::logERROR, (format("\nGWBSE ERR ...")).str());
                gwbse_logger.setPreface(ctp::logWARNING, (format("\nGWBSE WAR ...")).str());
                gwbse_logger.setPreface(ctp::logDEBUG, (format("\nGWBSE DBG ...")).str());

                // Initialize with options
                _gwbse.Initialize(&_gwbse_options);

                /* Only execute GWBSE if excited state is requested. This is a bit
                 * weird construction to have ground state calculation treated in
                 * exactly the same way for polarized QMMM.
                 */
                if (_state > 0) {
                    CTP_LOG(ctp::logDEBUG, *_log) << "Excited state via GWBSE: " << flush;
                    CTP_LOG(ctp::logDEBUG, *_log) << "  --- type:              " << _type << flush;
                    CTP_LOG(ctp::logDEBUG, *_log) << "  --- state:             " << _state << flush;
                    if (_has_overlap_filter) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: overlap  " <<  flush;
                    }
                    if (_has_osc_filter) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: osc.str. > " << _osc_threshold << flush;
                    }
                    if (_has_dQ_filter) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: crg.trs. > " << _dQ_threshold << flush;
                    }
                    if (_has_loc_filter){
                        if (_localiseonA){
                         CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: localisation on A > " << _loc_threshold << flush;
                        }else{
                            CTP_LOG(ctp::logDEBUG, *_log) << "  --- filter: localisation on B > " << _loc_threshold << flush;
                        }
                    }

                    if (_has_osc_filter && _has_dQ_filter) {
                        CTP_LOG(ctp::logDEBUG, *_log) << "  --- WARNING: filtering for optically active CT transition - might not make sense... " << flush;
                    }

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

                    // quasiparticle filter
                    if (_has_overlap_filter) {
                        if (iter == 0) {

                            // One - to - One LIST in 0th iteration
                            for (unsigned _i = 0; _i < orb_iter_input.QPdiagEnergies().size(); _i++) {
                                _state_index.push_back(_i);
                            }

                        } else {
                            // get AO overlap matrix
                            AOOverlap _dftoverlap;
                            // load dft  basis set (element-wise information) from xml file
                            BasisSet dftbs;
                            dftbs.LoadBasisSet(orb_iter_input.getDFTbasis());
                            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb_iter_input.getDFTbasis() << flush;

                            // fill auxiliary GW AO basis by going through all atoms
                            AOBasis dftbasis;
                            dftbasis.AOBasisFill(dftbs, orb_iter_input.QMAtoms());
                            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Basis of size " << dftbasis.AOBasisSize() << flush;

                            // Fill overlap
                            _dftoverlap.Fill(dftbasis);
                            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: " << _dftoverlap.Matrix().rows() << flush;


                            // 'LAMBDA' matrix of the present iteration
                            Eigen::MatrixXd lambda_N = orb_iter_input.LambdaMatrixQuasiParticle();

                            // 'LAMBDA' matrix of the previous iteration
                            string runFolder_N_1 = jobFolder + "/iter_" + boost::lexical_cast<string>(iter - 1);
                            string orbfile_N_1 = runFolder_N_1 + "/system.orb";
                            Orbitals _orbitals_N_1;
                            // load the QM data from serialized orbitals object
                           
                            CTP_LOG(ctp::logDEBUG, *_log) << " Loading QM data from " << orbfile_N_1 << flush;
                            _orbitals_N_1.ReadFromCpt(orbfile_N_1);
                            
                            Eigen::MatrixXd lambda_N_1 = _orbitals_N_1.LambdaMatrixQuasiParticle();
                            // calculate QP overlaps
                            
                            Eigen::MatrixXd qpoverlaps = lambda_N*_dftoverlap.Matrix()*lambda_N_1.transpose();

                            // test output
                            if (tools::globals::verbose) {
                                for (unsigned i = 0; i < qpoverlaps.rows(); i++) {
                                    for (unsigned j = 0; j < qpoverlaps.cols(); j++) {
                                        CTP_LOG(ctp::logDEBUG, *_log) << " [" << i << " , " << j << "]: " << qpoverlaps(i, j) << flush;
                                    }
                                }
                            }
                            
                            // filter for max absolute value (hopefully close to 1)
                            for (unsigned _j = 0; _j < qpoverlaps.cols(); _j++) {
                                int maxi = 0;
                                for (unsigned _i = 0; _i < qpoverlaps.rows(); _i++) {
                                    if (std::abs(qpoverlaps(_i, _j)) > std::abs(qpoverlaps(maxi, _j))) {
                                        maxi = _i;
                                    }
                                }
                                _state_index.push_back(maxi);
                                CTP_LOG(ctp::logDEBUG, *_log) << " [" << maxi << " , " << _j << "]: " << qpoverlaps(maxi, _j) << flush;
                            }

                        }

                    }

                    if (_has_osc_filter) {

                        // go through list of singlets
                        const std::vector<double>oscs = orb_iter_input.Oscillatorstrengths();
                        for (unsigned _i = 0; _i < oscs.size(); _i++) {

                            double osc = oscs[_i];
                            if (osc > _osc_threshold) _state_index.push_back(_i);
                        }

                    } else {
                        const VectorXfd & energies = (_type=="singlet") 
                        ? orb_iter_input.BSESingletEnergies() : orb_iter_input.BSETripletEnergies();

                        for (unsigned _i = 0; _i < energies.size(); _i++) {
                            _state_index.push_back(_i);
                        }
                    }


                    // filter according to charge transfer, go through list of excitations in _state_index
                    if (_has_dQ_filter) {
                        std::vector<int> _state_index_copy;
                        const std::vector< Eigen::VectorXd >& dQ_frag= (_type=="singlet") 
                        ? orb_iter_input.getFragmentChargesSingEXC():orb_iter_input.getFragmentChargesTripEXC();
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
                        ? orb_iter_input.getFragment_E_localisation_singlet():orb_iter_input.getFragment_E_localisation_triplet();
                        const std::vector< Eigen::VectorXd >& popH= (_type=="singlet") 
                        ? orb_iter_input.getFragment_H_localisation_singlet():orb_iter_input.getFragment_H_localisation_triplet();
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
                        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " WARNING: FILTER yielded no state. Taking lowest excitation" << flush;
                        _state_index.push_back(0);
                    }else{
                        if ( _type == "quasiparticle" ){
                            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded QP index: "<<_state_index[_state - 1 - orb_iter_input.getGWAmin()]<< flush;
                        }else {
                            CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Filter yielded state"<<_type<<":"<<_state_index[_state - 1]+1<< flush;
                        }
                    }
                    // - output its energy
                    if (_type == "singlet") {
                        energy___ex = orb_iter_input.BSESingletEnergies()[_state_index[_state - 1]] * tools::conv::hrt2ev; // to eV
                    } else if (_type == "triplet") {
                        energy___ex = orb_iter_input.BSETripletEnergies()[_state_index[_state - 1]] * tools::conv::hrt2ev; // to eV
                    } else if (_type == "quasiparticle") {
                        if ( _state > orb_iter_input.getNumberOfElectrons()  ) {
                            // unoccupied QPs: E_a = E_0 + eps_l
                            energy___ex = orb_iter_input.QPdiagEnergies()[_state_index[_state - 1 - orb_iter_input.getGWAmin()]] * tools::conv::hrt2ev; // to eV
                        } else {
                            // occupied QPs: E_c = E_0 - eps_h
                            energy___ex = -1.0*orb_iter_input.QPdiagEnergies()[_state_index[_state - 1 - orb_iter_input.getGWAmin()]] * tools::conv::hrt2ev; // to eV
                        }
                    }

                } // only if state >0

                if (!_static_qmmm) {
                    Density2Charges(_state_index);
                } // for polarized QMMM

            } //_do_gwbse


            /* new ESP fit only required for
             * - polarizable QMMM
             * AND
             * - GWBSE or DFT with internal DFTENGINE
             */
            if (!_static_qmmm && _qmpack->getPackageName() == "xtp" && !_do_gwbse) {
                Density2Charges();
            } // for polarized QMMM
            
            if(tools::globals::verbose){
              CTP_LOG(ctp::logDEBUG, *_log) <<"Calculated partial charges"<< flush;
              for(const QMAtom* atom:orb_iter_input.QMAtoms()){
                CTP_LOG(ctp::logDEBUG, *_log) <<atom->getType()<<" "<< atom->getPartialcharge()<< flush;
              }
            }

            // Test: go via GDMA instead of point charges, only for DFT with Gaussian!
            GDMA _gdma;
            if (_do_gdma) {
                if (_qmpack->getPackageName() != "gaussian" || _qmpack->getExecutable() != "g03") {

                    throw runtime_error(" Invalid QMPackage! " + _type + " Gaussian 03 only!");

                } else {
                    // get a GDMA object
                    _gdma.Initialize(&_gdma_options);
                    _gdma.setLog(_log);
                    _gdma.SetRunDir(runFolder);

                    CTP_LOG(ctp::logINFO, *_log) << "Running GDMA " << flush;
                    // prepare a GDMA input file
                    _gdma.WriteInputFile();

                    // run GDMA external
                    _gdma.RunExternal();

                    // parse output of gdma and update multipoles_full
                    _gdma.ParseOutputFile();

                } // use gdma
            } // _do_gdma


            assert(orb_iter_input.hasSelfEnergy());
            assert(orb_iter_input.hasQMEnergy());

            // EXTRACT & SAVE QM ENERGIES
            double energy___sf = orb_iter_input.getSelfEnergy();
            double energy_qmsf = orb_iter_input.getQMEnergy();
            double energy_qm__ = energy_qmsf - energy___sf;
            thisIter->setQMSF(energy_qm__, energy___sf, energy___ex);
            _job->setEnergy_QMMM(thisIter->getQMEnergy(), thisIter->getGWBSEEnergy(), thisIter->getSFEnergy(),
                    thisIter->getQMMMEnergy());

            // EXTRACT & SAVE QMATOM DATA
            std::vector< QMAtom* > &atoms = orb_iter_input.QMAtoms();

            thisIter->UpdatePosChrgFromQMAtoms(atoms, _job->getPolarTop()->QM0());

            if (_do_gdma) {

                // update PolarTop
                thisIter->UpdateMPSFromGDMA(_gdma.GetMultipoles(), _job->getPolarTop()->QM0());

            }

            // Update state variable
            if (_type == "quasiparticle" || _has_overlap_filter ){

                _state = _state_index[ _state -1 - orb_iter_input.getGWAmin() ] + 1 + orb_iter_input.getGWAmin();

            }



            // serialize this iteration
            if (_do_archive) {
                // save orbitals
                std::string ORB_FILE = runFolder + "/system.orb";
                CTP_LOG(ctp::logDEBUG, *_log) << "Archiving data to " << ORB_FILE << flush;
                orb_iter_input.WriteToCpt(ORB_FILE);
            }

            CTP_LOG(ctp::logINFO, *_log)
                    << format("Summary - iteration %1$d:") % (iterCnt + 1) << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... QM Size  = %1$d atoms") % int(atoms.size()) << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... E(QM)    = %1$+4.9e") % thisIter->getQMEnergy() << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... E(GWBSE) = %1$+4.9e") % thisIter->getGWBSEEnergy() << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... E(SF)    = %1$+4.9e") % thisIter->getSFEnergy() << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... E(FM)    = %1$+4.9e") % thisIter->getFMEnergy() << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... E(MM)    = %1$+4.9e") % thisIter->getMMEnergy() << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... E(QMMM)  = %1$+4.9e") % thisIter->getQMMMEnergy() << flush;
            if (!_static_qmmm) {
                CTP_LOG(ctp::logINFO, *_log)
                        << format("... RMS(dR)  = %1$+4.9e") % thisIter->getRMSdR() << flush;
                CTP_LOG(ctp::logINFO, *_log)
                        << format("... RMS(dQ)  = %1$+4.9e") % thisIter->getRMSdQ() << flush;
                CTP_LOG(ctp::logINFO, *_log)
                        << format("... SUM(dQ)  = %1$+4.9e") % thisIter->getSUMdQ() << flush;
            }
            // CLEAN DIRECTORY
            _qmpack->CleanUp();


            return 0;
        }


        template<class QMPackage>
        void QMMachine<QMPackage>::Density2Charges(std::vector<int> _state_index ){


                    // load DFT basis set (element-wise information) from xml file
                    BasisSet dftbs;
                    if (orb_iter_input.getDFTbasis() != "") {
                        dftbs.LoadBasisSet(orb_iter_input.getDFTbasis());
                        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Loaded DFT Basis Set " << orb_iter_input.getDFTbasis() << flush;
                    }

                    // fill DFT AO basis by going through all atoms
                    AOBasis dftbasis;
                    dftbasis.AOBasisFill(dftbs, orb_iter_input.QMAtoms());

                    Eigen::MatrixXd DMATGS = orb_iter_input.DensityMatrixGroundState();

                    Eigen::MatrixXd DMAT_tot = DMATGS; // Ground state + hole_contribution + electron contribution

                    if (_state > 0 ) {

                        if ( _type == "singlet" && _type == "triplet"){

                            std::vector<Eigen::MatrixXd > DMAT = orb_iter_input.DensityMatrixExcitedState(_type, _state_index[_state - 1]);
                            DMAT_tot = DMAT_tot - DMAT[0] + DMAT[1]; // Ground state + hole_contribution + electron contribution
                        } else if ( _type == "quasiparticle"){

                            Eigen::MatrixXd DMATQP = orb_iter_input.DensityMatrixQuasiParticle(  _state_index[_state - 1 - orb_iter_input.getGWAmin()]);

                            if ( _state > orb_iter_input.getNumberOfElectrons() ) {
                                DMAT_tot = DMAT_tot + DMATQP;
                            } else {
                                DMAT_tot = DMAT_tot - DMATQP;
                            }
                        }
                    }

                    // fill DFT AO basis by going through all atoms
                    std::vector< QMAtom* >& Atomlist = orb_iter_input.QMAtoms();

  Espfit esp = Espfit(_log);
                    

  int iter=_iters.size()-1;
  Eigen::MatrixXd DMAT_mixed;
  if (iter == 0) {
    DMAT_mixed = DMAT_tot;
  } else {
    DMAT_mixed = _alpha * _DMAT_old + (1 - _alpha) * DMAT_tot;
  }

  _DMAT_old = DMAT_mixed;
  esp.Fit2Density(Atomlist, DMAT_mixed, dftbasis,  "medium");

            return;
        }



        template<class QMPackage>
        QMMIter *QMMachine<QMPackage>::CreateNewIter() {

            QMMIter *newIter = new QMMIter(_iters.size());
            this->_iters.push_back(newIter);
            return newIter;
        }

        template<class QMPackage>
        bool QMMachine<QMPackage>::hasConverged() {

            _convg_dR = false;
            _convg_dQ = false;
            _convg_dE_QM = false;
            _convg_dE_MM = false;

            if (_iters.size() > 1) {

                QMMIter *iter_0 = _iters[_iters.size() - 2];
                QMMIter *iter_1 = _iters[_iters.size() - 1];

                double dR = iter_1->getRMSdR();
                double dQ = iter_1->getRMSdQ();
                double dE_QM = iter_1->getQMEnergy() - iter_0->getQMEnergy();
                double dE_MM = iter_1->getMMEnergy() - iter_0->getMMEnergy();

                CTP_LOG(ctp::logINFO, *_log)
                        << format("... dE_QM  = %1$+4.9e") % dE_QM << flush;
                CTP_LOG(ctp::logINFO, *_log)
                        << format("... dE_MM  = %1$+4.9e") % dE_MM << flush;

                if (dR <= _crit_dR) _convg_dR = true;
                if (dQ <= _crit_dQ) _convg_dQ = true;
                if (dE_QM * dE_QM <= _crit_dE_QM * _crit_dE_QM) _convg_dE_QM = true;
                if (dE_MM * dE_MM <= _crit_dE_MM * _crit_dE_MM) _convg_dE_MM = true;
            }

            _isConverged = ((_convg_dR && _convg_dQ) && (_convg_dE_QM && _convg_dE_MM));



            CTP_LOG(ctp::logINFO, *_log)
                    << format("... Convg dR = %s") % (_convg_dR ? "true" : "false") << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... Convg dQ = %s") % (_convg_dQ ? "true" : "false") << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... Convg QM = %s") % (_convg_dE_QM ? "true" : "false") << flush;
            CTP_LOG(ctp::logINFO, *_log)
                    << format("... Convg MM = %s") % (_convg_dE_MM ? "true" : "false") << flush;

            return _isConverged;
        }



        // REGISTER QM PACKAGES
        template class QMMachine<QMPackage>;



    }
}
