/* 
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "exciton.h"
//#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/gwbseengine.h>


using namespace std;

namespace votca {
    namespace xtp {

    


void Exciton::Initialize(Property* options) {

            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults( options, "xtp" );



            _do_optimize = false;
            _restart_opt = false; // GEOOPT Option
            // _do_guess=false; //Writing guess for dimer calculation
            
            string key = "options." + Identify();
            _archive_file = options->get(key + ".archive").as<string>();
            _reporting =  options->get(key + ".reporting").as<string> ();
            
            // job tasks
            string _tasks_string = options->get(key + ".tasks").as<string> ();
            if (_tasks_string.find("optimize") != std::string::npos) _do_optimize = true;
            
            // GWBSEENGINE options
            _gwbseengine_options = options->get(key + ".gwbse_engine"); // not sure if this works

            // XML Output needed at all?
            _xml_output = options->get(key + ".output").as<string> ();

            // options for dft package
            string _package_xml = options->get(key + ".dftpackage").as<string> ();
            load_property_from_xml(_package_options, _package_xml.c_str());
            key = "package";
            _package = _package_options.get(key + ".name").as<string> ();

            // MOLECULE properties
            key = "options." + Identify();
            _xyzfile = options->get(key + ".molecule").as<string>();
            
            // if optimization is chosen
            if ( _do_optimize ){
                  _opt_state = options->get(key + ".optimize.state").as<int> ();
                  _displacement = options->get(key + ".optimize.displacement").as<double> ();
                  _convergence = options->get(key + ".optimize.convergence").as<double>();
                  _spintype = options->get(key + ".optimize.spintype").as<string> ();
                  _trust_radius = options->get(key + ".optimize.trust").as<double> ();
                  _forces = options->get(key + ".optimize.forces").as<string> ();
                  _opt_type = options->get(key + ".optimize.type").as<string> ();
                  string _restart_string = options->get(key + ".optimize.restart").as<string> ();
                  if ( _restart_string == "restart" ) _restart_opt = true;
                  if ( _restart_opt  && ! boost::filesystem::exists("restart.opt")) {
                      cout << " ERROR: Restart requested but restart.opt file not found!" << endl;
                      exit(1);
                  };
            }

            // register all QM packages (Gaussian, TURBOMOLE, etc)
            QMPackageFactory::RegisterAll();


           
            
}

bool Exciton::Evaluate() {

    
    if ( _reporting == "silent")   _log.setReportLevel( ctp::logERROR ); // only output ERRORS, GEOOPT info, and excited state info for trial geometry
    if ( _reporting == "noisy")    _log.setReportLevel( ctp::logDEBUG ); // OUTPUT ALL THE THINGS
    if ( _reporting == "default")  _log.setReportLevel( ctp::logINFO );  // 
    
    _log.setMultithreading( true );
    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ..."); 

    ctp::TLogLevel _ReportLevel = _log.getReportLevel( ); // backup report level
    
    if ( _do_optimize ){
        _trust_radius = _trust_radius * tools::conv::ang2bohr; // initial trust radius in a.u.
        _trust_radius_max = 0.1; 
        _displacement = _displacement * tools::conv::ang2bohr; // initial trust radius in a.u.
        _log.setReportLevel( ctp::logINFO );
        CTP_LOG(ctp::logINFO,_log) << "Requested geometry optimization of excited state " << _spintype << " " << _opt_state << flush; 
        CTP_LOG(ctp::logINFO,_log) << "... forces calculation using " << _forces << " differences " << flush;
        CTP_LOG(ctp::logINFO,_log) << "... displacement for numerical forces:           " << _displacement << " Bohr" << flush; 
        CTP_LOG(ctp::logINFO,_log) << "... convergence of total energy:                 " << _convergence << " Hartree" << flush;
        CTP_LOG(ctp::logINFO,_log) << "... initial trust radius:                        " << _trust_radius << " Bohr " << flush;
        _log.setReportLevel( _ReportLevel );
    }
    
    // Get orbitals object
    Orbitals _orbitals;
    
    // Read molecular geometry from xyz file and store in a segment (WHY SEGMENT?)
    std::vector <ctp::Segment* > _segments;
    ctp::Segment _segment(0, "mol");
    CTP_LOG(ctp::logDEBUG,_log) << " Reading molecular coordinates from " << _xyzfile << flush;
    ReadXYZ( &_segment, _xyzfile );
    _segments.push_back(&_segment);
    
    // Get and initialize QMPackage for DFT ground state
    QMPackage *_qmpackage =  QMPackages().Create( _package );
    _qmpackage->setLog( &_log );       
    _qmpackage->Initialize( &_package_options );
    _qmpackage->setRunDir(".");

    // Get GWBSEENGINE Object and initialize
    GWBSEENGINE _gwbse_engine;
    _gwbse_engine.setLog( &_log );
    _gwbse_engine.Initialize(&_gwbseengine_options, _archive_file);
    
    // Run GWBSE
    _gwbse_engine.ExcitationEnergies( _qmpackage, _segments, &_orbitals);





    /*
    //_iteration = 0;
    string _geoopt_preface = "\n... ... GEOOPT " + boost::lexical_cast<std::string>(_iteration);
    if ( _do_optimize ) _log.setPreface(ctp::logINFO, _geoopt_preface );
    
    if ( ! _restart_opt ){
    
    // do GWBSE for start geometry
    if ( _reporting == "silent" ) _log.setReportLevel( ctp::logINFO );
    
    _gwbse_engine.ExcitationEnergies( _qmpackage, _segments, &_orbitals);

    }

   // only do this, if optimization is requested 
   if ( _do_optimize ){
   
       GeometryOptimization _geoopt;
       
       
       
      // optimization implies running DFT for displaced geometries
      _do_dft_input = true;
      _do_dft_run = true;
      _do_dft_parse = true;
      
      // convergence checks
      bool _converged          = false;
      double energy=0;
      if ( _restart_opt ){
          // reload data from restart.opt
          ReloadState();
          // exit(0);
      } else {

         // get total energy for this excited state
          energy = _orbitals.GetTotalEnergy( _spintype, _opt_state );
         CTP_LOG(ctp::logINFO,_log) << " Total energy at initial geometry for " << _spintype << " " << _opt_state << " is " << energy << " Hartree " << flush;
         CTP_LOG(ctp::logINFO,_log) << "   " << flush;
         _log.setReportLevel( _ReportLevel ); // go silent again during numerical force calculation, if requested
      
      }
      //double energy_new;
      // now iterate
      while ( ! _converged ){ 
           _geoopt_preface = "\n... ... GEOOPT " + boost::lexical_cast<std::string>(_iteration+1);
           _log.setPreface(ctp::logINFO, _geoopt_preface );
           
           
           std::vector <ctp::Segment* > _molecule;
           ctp::Segment _current_coordinates(0, "mol");
           if ( _restart_opt) {
               // fill current coordinates from reloaded data
               CTP_LOG(ctp::logDEBUG,_log) << " Reading molecular coordinates from " << _xyzfile << flush;
               ReadXYZ( &_current_coordinates, _xyzfile ); // to initialize
               Coord2Segment( &_current_coordinates); // overwrite coordinates
               _molecule.push_back( &_current_coordinates );
                ExcitationEnergies( _qmpackage, _molecule, &_orbitals);
                energy =  _orbitals.GetTotalEnergy(_spintype, _opt_state);
               _restart_opt = false;
           } else {
               // get current coordinates from Orbitals object
               Orbitals2Segment( &_current_coordinates , &_orbitals);
               _molecule.push_back( &_current_coordinates );
           }
          
           
           // displace all atoms in each cartesian coordinate and get new energy
           // for each atom
           std::vector< ctp::Atom* > _atoms;
           std::vector< ctp::Atom* > ::iterator ait;
           _atoms = _current_coordinates.Atoms();
           
           // initialize the Geometry Optimization objects
           if ( _geoopt.Iteration() == 0 ) _geoopt.Initialize(_atoms.size1()); 
           
           // calculate numerical forces
           if ( _forces == "forward" ) NumForceForward( energy, _atoms, _force, _qmpackage, _molecule, &_orbitals );
           if ( _forces == "central" ) NumForceCentral( energy, _atoms, _force, _qmpackage, _molecule, &_orbitals );

           _geoopt.Checkpoint(_molecule); 
           
           if ( _opt_type != "md" ){
           
           _step_accepted = false; // at every iteration
           
           // initial trial -> update Hessian
           _update_hessian = true;
           double _RMSForce;
           double _MaxForce;
           double _RMSStep;
           double _MaxStep;
           double energy_delta;
           double energy_new;
           double _old_TR = _trust_radius;
           while ( ! _step_accepted ){
           
               _geoopt.BFGSStep(_iteration, _update_hessian, _force,_force_old, _current_xyz, _old_xyz, _hessian,  _xyz_shift, _trial_xyz);
           
               // put trial coordinates into _segment
               int _i_atom = 0;
               for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                   // put trial coordinates (_trial_xyz is in Bohr, segments in nm)
                   vec _pos_displaced( _trial_xyz(_i_atom,0)*tools::conv::bohr2nm , _trial_xyz(_i_atom,1)*tools::conv::bohr2nm, _trial_xyz(_i_atom,2)*tools::conv::bohr2nm );
                   (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment
                   _i_atom++;
               }

               // for the updated geometry, get new reference energy
               if ( _reporting == "silent" ) _log.setReportLevel( ctp::logINFO );
               ExcitationEnergies(_qmpackage, _molecule, &_orbitals);
               energy_new  = _orbitals.GetTotalEnergy(_spintype, _opt_state );
               
               // CTP_LOG(ctp::logINFO,_log) << " Total energy of " << _spintype << " " << _opt_state << " is " << energy_new << " a.u. at " << _iteration << flush;
               energy_delta = energy_new - energy;
               _old_TR = _trust_radius;
               // check the energy at the trial coordinates
               if ( energy_delta > 0.0 ){
                   // total energy has unexpectedly increased, half the trust radius
                   _trust_radius = 0.5 * _trust_radius;
                   // Hessian should not be updated for reduced trust radius
                   _update_hessian = false;
                   CTP_LOG(ctp::logINFO,_log) << " BFGS-TRM: total energy increased, this step is rejected and the new trust radius is: " << _trust_radius << endl;
                   // repeat BGFSStep with new trust radius
               } else {
                   // total energy has decreased, we accept the step but might update the trust radius
                   _step_accepted = true;
                   CTP_LOG(ctp::logINFO,_log) << " BFGS-TRM: step accepted" << flush;
                   // adjust trust radius, if required
                   double _tr_check = energy_delta / _delta_energy_estimate;
                   if ( _tr_check > 0.75 && 1.25*sqrt(_norm_delta_pos) > _trust_radius ){
                      _trust_radius = 2.0 * _trust_radius;
                      // if ( _trust_radius > _trust_radius_max) _trust_radius = _trust_radius_max;
		      CTP_LOG(ctp::logINFO,_log) << " BFGS-TRM: increasing trust radius to " << _trust_radius << endl; 
                   } else if ( _tr_check < 0.25 ) {
                      _trust_radius = 0.25*sqrt(_norm_delta_pos);
                      CTP_LOG(ctp::logINFO,_log) << " BFGS-TRM: reducing trust radius to " << _trust_radius << endl;
                   }
               
		   }
           }
           
           // after the step is accepted, we can shift the stored data   
           energy       = energy_new;
           _old_xyz     = _current_xyz;
           _current_xyz = _trial_xyz;
           _force_old   = _force; 

           // checking convergence
           bool _energy_converged   = false;
           bool _RMSForce_converged = false;
           bool _MaxForce_converged = false;
           bool _RMSStep_converged  = false;
           bool _MaxStep_converged  = false;

           _xyz_shift = _current_xyz - _old_xyz;
           
           _RMSForce = linalg_getRMS( _force );
           _MaxForce = linalg_getMax( _force,true );
           _RMSStep  = linalg_getRMS( _xyz_shift );
           _MaxStep  = linalg_getMax( _xyz_shift,true );


           
           if ( std::abs(energy_delta) < _convergence) _energy_converged = true;
           if ( std::abs(_RMSForce) < 0.01 ) _RMSForce_converged = true;
           if ( std::abs(_MaxForce) < 0.01 ) _MaxForce_converged = true;
           if ( std::abs(_RMSStep ) < 0.01 ) _RMSStep_converged = true;
           if ( std::abs(_MaxStep) < 0.01 ) _MaxStep_converged = true;
           
           if ( _energy_converged && _RMSForce_converged && _MaxForce_converged && _RMSStep_converged && _MaxStep_converged ) _converged = true;

        // dump information to restart file
        
        FILE *restart;
        string FILE = "restart.opt";
        restart = fopen(FILE.c_str(),"w");
        fprintf(restart, "atoms %d \n", _natoms);
        fprintf(restart, "iteration %d \n", _iteration);
        // output current xyz
        fprintf(restart, "current xyz \n");
        for ( int _i_atom=0; _i_atom< _natoms; _i_atom++){
            fprintf(restart, "%le %le %le \n", _current_xyz(_i_atom,0) , _current_xyz(_i_atom,1) , _current_xyz(_i_atom,2)  );
        }
        // output old xyz
        fprintf(restart, "old xyz \n");
        for ( int _i_atom=0; _i_atom< _natoms; _i_atom++){
            fprintf(restart, "%le %le %le \n", _old_xyz(_i_atom,0) , _old_xyz(_i_atom,1) , _old_xyz(_i_atom,2)  );
        }
        // output current force
        fprintf(restart, "current force \n");
        for ( int _i_atom=0; _i_atom< _natoms; _i_atom++){
            fprintf(restart, "%le %le %le \n", _force(_i_atom,0) , _force(_i_atom,1) , _force(_i_atom,2)  );
        }
        // output old force
        fprintf(restart, "old force \n");
        for ( int _i_atom=0; _i_atom< _natoms; _i_atom++){
            fprintf(restart, "%le %le %le \n", _force_old(_i_atom,0) , _force_old(_i_atom,1) , _force_old(_i_atom,2)  );
        }
        // output Hessian
        fprintf(restart, "Hessian \n");
        for ( unsigned _i = 0; _i < _hessian.size1() ; _i++){
            for ( unsigned _j = 0; _j < _hessian.size2() ; _j++){
                
                fprintf(restart, "%d %d %le \n", _i,_j,_hessian(_i,_j)); 
                
            }
        }
 
        fclose(restart);

           CTP_LOG(ctp::logINFO,_log) << " Summary of iteration " << _iteration +1 << flush;
           CTP_LOG(ctp::logINFO,_log) << "    total energy of " << _spintype << " " << _opt_state << " is " << energy_new << " Hartree "  << flush;
           CTP_LOG(ctp::logINFO,_log) << "    BFGS-TRM: trust radius was " << _old_TR << " Bohr" << flush;
           CTP_LOG(ctp::logINFO,_log) << "    convergence: " << flush;
           CTP_LOG(ctp::logINFO,_log) << "       energy change: " << setprecision(12) << energy_delta << " a.u. " << Convergence( _energy_converged ) << flush;
           CTP_LOG(ctp::logINFO,_log) << "       RMS force:     " << setprecision(12) << _RMSForce << Convergence( _RMSForce_converged ) << flush;
           CTP_LOG(ctp::logINFO,_log) << "       Max force:     " << setprecision(12) << _MaxForce << Convergence( _MaxForce_converged ) << flush;
           CTP_LOG(ctp::logINFO,_log) << "       RMS step:      " << setprecision(12) << _RMSStep << Convergence( _RMSStep_converged ) << flush;
           CTP_LOG(ctp::logINFO,_log) << "       Max step:      " << setprecision(12) << _MaxStep << Convergence( _MaxStep_converged ) << flush;
           CTP_LOG(ctp::logINFO,_log) << "       Geometry       " << Convergence( _converged ) << flush;
           CTP_LOG(ctp::logINFO,_log) << "   " << flush;
           _log.setReportLevel( _ReportLevel ); // go silent again
           
           
           _iteration++;
	   // exit(0);
	   } else {
             _converged = true;
	   }
       } // optimization while loop
       
   }// if optimization
  */
    
    
    
    
       CTP_LOG(ctp::logDEBUG,_log) << "Saving data to " << _archive_file << flush;
       std::ofstream ofs( ( _archive_file).c_str() );
       boost::archive::binary_oarchive oa( ofs );


       
       oa << _orbitals;
       ofs.close();

    //Property *_job_output = &_summary.add("output","");
    //tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");
    // CTP_LOG(ctp::logDEBUG,_log) << "Writing output to " << _xml_output << flush;
    //std::ofstream ofout (_xml_output.c_str(), std::ofstream::out);
    //if ( _do_gwbse ){
    //    ofout << (_summary.get("output"));    
   // }
    //ofout.close();
    
    return true;
}





/* void Exciton::Coord2Segment( ctp::Segment* _segment){
    
    
            std::vector< ctp::Atom* > _atoms;
            std::vector<  ctp::Atom* > ::iterator ait;
            _atoms = _segment->Atoms();
            
            string type;
            int _i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                vec position(_current_xyz(_i_atom,0)*tools::conv::bohr2nm, _current_xyz(_i_atom,1)*tools::conv::bohr2nm, _current_xyz(_i_atom,2)*tools::conv::bohr2nm); // current_xyz has Bohr, votca stores nm
                (*ait)->setQMPos(position);
                _i_atom++;
               
            }
            
            return;
} */


/* reload previous optimization state*/
/* void Exciton::ReloadState(){
    
    // open restart.opt
    ifstream in;
    string _filename = "restart.opt";
    CTP_LOG(ctp::logDEBUG,_log) << " Reloading from restart.opt " << flush;
    in.open(_filename.c_str(), ios::in);
    if (!in) throw runtime_error(string("Error reading coordinates from restart.opt"));
    string label;
    string type;
    // atoms
    in >> label;
    if ( label != "atoms" ) throw runtime_error(string("Expected atoms in restart.opt"));
    in >> _natoms;
    // initialize matrices
    _force       = ub::zero_matrix<double>(_natoms ,3);
    _force_old   = ub::zero_matrix<double>(_natoms ,3);
    _xyz_shift   = ub::zero_matrix<double>(_natoms ,3);
    _current_xyz = ub::zero_matrix<double>(_natoms ,3);
    _old_xyz     = ub::zero_matrix<double>(_natoms ,3);
    _trial_xyz   = ub::zero_matrix<double>(_natoms ,3);
    _hessian     = ub::zero_matrix<double>(3*_natoms,3*_natoms);
    
    // info about last iteration
    in >> label;
    if ( label != "iteration" ) throw runtime_error(string("Expected iteration in restart.opt"));
    in >> _iteration;
    
    // current xyz
    in >> label; 
    if ( label != "current" ) throw runtime_error(string("Expected current xyz in restart.opt"));
    in >> type;
    if ( type != "xyz" ) throw runtime_error(string("Expected current xyz in restart.opt"));
    int _i_atom = 0;
    while ( _i_atom < _natoms ){
        in >> _current_xyz(_i_atom,0);
        in >> _current_xyz(_i_atom,1);
        in >> _current_xyz(_i_atom,2);
        _i_atom++;
    }

    // old xyz
    in >> label; 
    if ( label != "old" ) throw runtime_error(string("Expected old xyz in restart.opt"));
    in >> type;
    if ( type != "xyz" ) throw runtime_error(string("Expected old xyz in restart.opt"));
    _i_atom = 0;
    while ( _i_atom < _natoms ){
        in >> _old_xyz(_i_atom,0);
        in >> _old_xyz(_i_atom,1);
        in >> _old_xyz(_i_atom,2);
        _i_atom++;
    }
    
    // current force
    in >> label; 
    if ( label != "current" ) throw runtime_error(string("Expected current force in restart.opt"));
    in >> type;
    if ( type != "force" ) throw runtime_error(string("Expected current force in restart.opt"));
    _i_atom = 0;
    while ( _i_atom < _natoms ){
        in >> _force(_i_atom,0);
        in >> _force(_i_atom,1);
        in >> _force(_i_atom,2);
        _i_atom++;
    }   
    
    // old force
    in >> label; 
    if ( label != "old" ) throw runtime_error(string("Expected old force in restart.opt"));
    in >> type;
    if ( type != "force" ) throw runtime_error(string("Expected old force in restart.opt"));
    _i_atom = 0;
    while ( _i_atom < _natoms ){
        in >> _force_old(_i_atom,0);
        in >> _force_old(_i_atom,1);
        in >> _force_old(_i_atom,2);
        _i_atom++;
    }   
       
    // hessian
    int _i_in;
    int _j_in;
    in >> label; 
    if ( label != "Hessian" ) throw runtime_error(string("Expected Hessian in restart.opt"));
    for ( int _i =0 ; _i < 3*_natoms; _i++){
        for ( int _j =0 ; _j < 3*_natoms; _j++){
            in >> _i_in;
            in >> _j_in;
            in >> _hessian(_i,_j);
        }
    }
    
    return;
} */





/* Calculate forces on atoms numerically by central differences */
/*void Exciton::NumForceCentral(double energy, std::vector<ctp::Atom*> _atoms, ub::matrix<double>& _force, 
        QMPackage* _qmpackage, std::vector<ctp::Segment*> _molecule, Orbitals* _orbitals){
    

    std::vector< ctp::Atom* > ::iterator ait;
    int _i_atom = 0;
    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
       // get its current coordinates
       vec _current_pos = (*ait)->getQMPos(); // in nm
    
       // go through all cartesian components
       for ( int _i_cart = 0 ; _i_cart < 3; _i_cart++ ){
           
           // get displacement vector in positive direction
           vec _displaced(0, 0, 0);
           if (_i_cart == 0) {
              _displaced.setX(_displacement / tools::conv::nm2bohr ); // x, _displacement in Bohr
              _current_xyz(_i_atom,_i_cart ) = _current_pos.getX() * tools::conv::nm2bohr; // in Bohr
           }
           if (_i_cart == 1) {
              _current_xyz(_i_atom,_i_cart ) = _current_pos.getY() * tools::conv::nm2bohr; // in Bohr
              _displaced.setY(_displacement / tools::conv::nm2bohr ); // y, _displacement in in Angstrom
           }
           if (_i_cart == 2) {
              _current_xyz(_i_atom,_i_cart ) = _current_pos.getZ() * tools::conv::nm2bohr; // in Bohr
              _displaced.setZ(_displacement / tools::conv::nm2bohr ); // z, _displacement in in Angstrom
           }
                            
           // update the coordinate
           vec _pos_displaced = _current_pos + _displaced;
           (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

           // run DFT and GW-BSE for this geometry
           ExcitationEnergies(_qmpackage, _molecule, _orbitals);
                            
           // get total energy for this excited state
           double energy_displaced_plus = _orbitals->GetTotalEnergy(_spintype, _opt_state);
           
           // get displacement vector in negative direction

           // update the coordinate
           _pos_displaced = _current_pos - 2.0 * _displaced;
           (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

           // run DFT and GW-BSE for this geometry
           ExcitationEnergies(_qmpackage, _molecule, _orbitals);
                            
           // get total energy for this excited state
           double energy_displaced_minus =  _orbitals->GetTotalEnergy(_spintype, _opt_state);
           
           // calculate force and put into matrix
           _force(_i_atom,_i_cart) = 0.5 * (energy_displaced_minus - energy_displaced_plus) / _displacement ; // force a.u./a.u.
                                                  
           (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
        }
               
      _i_atom++;
   }
    

    
    return;
} */



/* Calculate forces on atoms numerically by forward differences */
/* void Exciton::NumForceForward(double energy, std::vector<ctp::Atom*> _atoms, ub::matrix<double>& _force, 
        QMPackage* _qmpackage, std::vector<ctp::Segment*> _molecule, Orbitals* _orbitals){

    
    // check if file "forces.resume" exists
    string _filename = "forces.resume";
    bool _forcefile_exists = boost::filesystem::exists(_filename);
    ifstream in;
    bool _forces_resume = false;
    double f_in;
    if ( _forcefile_exists ) {
        // get iteration number
        in.open(_filename.c_str(), ios::in);
        int _iter_resume;
        in >> _iter_resume;
        if ( _iter_resume != _iteration ){
            CTP_LOG(ctp::logDEBUG,_log) << " Forces in " << _filename << " from iteration " << _iter_resume << " and not iteration " << _iteration << flush;
            CTP_LOG(ctp::logDEBUG,_log) << "  discard data..." << flush;
        } else {
            CTP_LOG(ctp::logDEBUG,_log) << " Resuming force calculation from  " << _filename << flush;
            _forces_resume = true;
            in >> f_in;
        }
    }



    ofstream out;
    vector< ctp::Atom* > ::iterator ait;
    int _i_atom = 0;
    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
        

          // get its current coordinates
          vec _current_pos = (*ait)->getQMPos(); // in nm
          //       if ( _i_atom < 2 ) {        
          for ( int _i_cart = 0 ; _i_cart < 3; _i_cart++ ){

             // try reading from forces.resume
             if ( _forces_resume && !in.eof()){
                 cout << " reading " << _i_atom << " " << _i_cart << endl;
                _force(_i_atom,_i_cart) = f_in;
                in >> f_in;
                if (_i_cart == 0) {
                   _current_xyz(_i_atom,_i_cart ) = _current_pos.getX() *  tools::conv::nm2bohr; // in Bohr
                }
                if (_i_cart == 1) {
                 _current_xyz(_i_atom,_i_cart ) = _current_pos.getY() *  tools::conv::nm2bohr; // in Bohr
                }
                if (_i_cart == 2) {
                 _current_xyz(_i_atom,_i_cart ) = _current_pos.getZ() *  tools::conv::nm2bohr; // in Bohr
                }
             } else {
        
                 if ( in.is_open() ) {
                     // close for reading
                     in.close();
                     out.open( _filename.c_str(), ios::app  );
                 } else {
                     // and open to append
                     if ( ! out.is_open() ) {
                         out.open( _filename.c_str(), ios::app  );
                         out << _iteration << endl;
                     }
                 }    
                 cout << " calculating " << _i_atom << " " << _i_cart << endl;
              
              // get displacement std::vector
              vec _displaced(0, 0, 0);
              if (_i_cart == 0) {
                 _displaced.setX(_displacement / tools::conv::nm2bohr ); // x, _displacement in Bohr
                 _current_xyz(_i_atom,_i_cart ) = _current_pos.getX()*tools::conv::nm2bohr; // in Bohr
              }
              if (_i_cart == 1) {
                 _current_xyz(_i_atom,_i_cart ) = _current_pos.getY()* tools::conv::nm2bohr; // in Bohr
                 _displaced.setY(_displacement /  tools::conv::nm2bohr ); // y, _displacement in in Angstrom
              }
              if (_i_cart == 2) {
                 _current_xyz(_i_atom,_i_cart ) = _current_pos.getZ() *  tools::conv::nm2bohr; // in Bohr
                 _displaced.setZ(_displacement /  tools::conv::nm2bohr ); // z, _displacement in in Angstrom
              }
                            
              // update the coordinate
              vec _pos_displaced = _current_pos + _displaced;
              (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

              // run DFT and GW-BSE for this geometry
              ExcitationEnergies(_qmpackage, _molecule, _orbitals);
                            
              // get total energy for this excited state
              double energy_displaced =  _orbitals->GetTotalEnergy(_spintype, _opt_state);
              // calculate force and put into matrix
           
              _force(_i_atom,_i_cart) = (energy - energy_displaced) / _displacement ; // force a.u./a.u.
              cout << "writing force " << _i_atom << " " << _i_cart << endl;
              out << setprecision(12) << _force(_i_atom,_i_cart) << endl;
              
              
              
              (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
           }
          }  
      _i_atom++;
   
    }
    
    if ( out.is_open() ) out.close();
    return;
} */


void Exciton::ReadXYZ(ctp::Segment* _segment, string filename){
    
                string line;
                std::ifstream in;
                
           
                string label, type;
                vec pos;

                
                in.open(filename.c_str(), std::ios::in);
                if (!in) throw runtime_error(string("Error reading coordinates from: ")
                        + filename);

                
                int atomCount = 1;

                if (in.is_open() ) {
                    while ( in.good() ) {
                    std::getline(in, line);

                    std::vector< string > split;
                    Tokenizer toker(line, " \t");
                    toker.ToVector(split);
                    if ( !split.size()      ||
                          split.size() != 4 ||
                          split[0] == "#"   ||
                          split[0].substr(0,1) == "#" ) { continue; }

            // Interesting information written here: e.g. 'C 0.000 0.000 0.000'
                        atomCount++;
                        string element = split[0];
                        double x = boost::lexical_cast<double>( split[1] ) / 10.; //°A to NM
                        double y = boost::lexical_cast<double>( split[2] ) / 10.;
                        double z = boost::lexical_cast<double>( split[3] ) / 10.;
                        vec Pos = vec(x,y,z);
                        ctp::Atom *pAtom = new ctp::Atom(atomCount, element);
                        pAtom->setPos(Pos);
                        pAtom->setQMPart(atomCount, Pos);
                        pAtom->setElement(element);
                        _segment->AddAtom(pAtom);
    
        }
    }
    else {
        throw std::runtime_error("No such file: '"+filename+"'.");
    }
          
return;    
}

/* void Exciton::Orbitals2Segment(ctp::Segment* _segment, Orbitals* _orbitals){
    
            std::vector< ctp::QMAtom* > _atoms;
            std::vector< ctp::QMAtom* > ::iterator ait;
            _atoms = _orbitals->QMAtoms();
            
            string type;
            int id = 1;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                
                // Atom *pAtom = new Atom(id++, type);
                
                type = (*ait)->type;
                double x = (*ait)->x;
                double y = (*ait)->y;
                double z = (*ait)->z;
                ctp::Atom *pAtom = new ctp::Atom(id++, type);
                // cout << type << " " << x << " " << y << " " << z << endl;
                vec position(x/10, y/10, z/10); // xyz has Angstrom, votca stores nm
                pAtom->setPos(position);
                pAtom->setQMPart(id, position);
                pAtom->setElement(type);
                _segment->AddAtom(pAtom);
            }
            return;
} */




  
    
    }
}
