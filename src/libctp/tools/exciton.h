/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef _VOTCA_CTP_EXCITON_H
#define _VOTCA_CTP_EXCITON_H

#include <stdio.h>

#include <votca/ctp/logger.h>
#include <votca/ctp/gwbse.h>
#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/segment.h>

namespace votca { namespace ctp {
    using namespace std;
    
class Exciton : public QMTool
{
public:

    Exciton() { };
   ~Exciton() { };

    string Identify() { return "exciton"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    GWBSE _gwbse;
 

private:
    
    string      _orbfile;
    string      _logfile;
    string      _xyzfile;

    string      _package;
    Property    _package_options;
    Property    _gwbse_options;
    
    string      _output_file;
    
    Logger      _log;
    
    bool _do_dft_input;
    bool _do_dft_run;
    bool _do_dft_parse;
    bool _do_gwbse;
    bool _do_optimize;
    
    int _opt_state;
    double _displacement;
    double _convergence;
    string _spintype;
    
    
    int _natoms;
    
    ub::matrix<double> _force;
    ub::matrix<double> _force_old;
    ub::matrix<double> _xyz_shift;
    ub::matrix<double> _speed;
    ub::matrix<double> _current_xyz;
    ub::matrix<double> _old_xyz; 
    
    void ExcitationEnergies( QMPackage* _qmpackage, vector <Segment* > _segments, Orbitals* _orbitals );
    void ReadXYZ( Segment* _segment, string filename);
    void Orbitals2Segment(Segment* _segment, Orbitals* _orbitals);
    double GetTotalEnergy( Orbitals* _orbitals, string _spintype, int _opt_state );
    
    
    void ConjugateGradientStep( int& _iteration, ub::matrix<double>& _force, ub::matrix<double>& _force_old,  ub::matrix<double>& _xyz_shift, ub::matrix<double>& _speed, ub::matrix<double>& _current_xyz, ub::matrix<double>&  _old_xyz );

};

void Exciton::Initialize(Property* options) {

            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults(options);

            _do_dft_input = false;
            _do_dft_run = false;
            _do_dft_parse = false;
            _do_gwbse = false;
            _do_optimize = false;
            
            string key = "options." + Identify();
            _output_file = options->get(key + ".archive").as<string>();

            // options for GWBSE package
            string _gwbse_xml = options->get(key + ".gwbse").as<string> ();
            //cout << endl << "... ... Parsing " << _package_xml << endl ;
            load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());

            
            // job tasks
            string _tasks_string = options->get(key + ".tasks").as<string> ();
            if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
            if (_tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
            if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;
            if (_tasks_string.find("optimize") != std::string::npos) _do_optimize = true;
            
            // something unique
            // _package = options->get(key + ".package").as<string> ();

            // options for dft package
            string _package_xml = options->get(key + ".package").as<string> ();
            //cout << endl << "... ... Parsing " << _package_xml << endl ;
            load_property_from_xml(_package_options, _package_xml.c_str());
            key = "package";
            _package = _package_options.get(key + ".name").as<string> ();
            
            
            key = "options." + Identify();
            // orbitals file or pure DFT output
            _orbfile = options->get(key + ".molecule.orbitals").as<string> ();
            _logfile = options->get(key + ".molecule.log").as<string> ();

            // test for parsing an xyz file
            if ( _do_dft_input ) {
                _xyzfile = options->get(key + ".molecule.xyz").as<string>();
            }
            
            // if optimization is chosen
            if ( _do_optimize ){
                  _opt_state = options->get(key + ".optimize.state").as<int> ();
                  _displacement = options->get(key + ".optimize.displacement").as<double> ();
                  _convergence = options->get(key + ".optimize.convergence").as<double>();
                  _spintype = options->get(key + ".optimize.spintype").as<string> ();
            }

            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
            string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/qmpackages/") + _package + string("_idft_pair.xml");
            //load_property_from_xml(_package_options, xmlFile);

            // register all QM packages (Gaussian, TURBOMOLE, etc)
            QMPackageFactory::RegisterAll();
            
}

bool Exciton::Evaluate() {

    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    
    
    if ( _do_optimize ){
       LOG(logDEBUG,_log) << "Requested geometry optimization of excited state " << _spintype << " " << _opt_state << flush; 
       LOG(logDEBUG,_log) << " ... displacement for numerical forces:          " << _displacement << " Angstrom" << flush; 
       LOG(logDEBUG,_log) << " ... convergence of total energy (?):            " << _convergence << " Ryd. (?)" << flush; 
       LOG(logDEBUG,_log) << " ... GOOD LUCK YOU CRAZY BASTARD!            "  << flush; 
    }
    
    
    vector <Segment* > _segments;
    // Create a new segment
    Segment _segment(0, "mol");
                
    if (_do_dft_input) {
       ReadXYZ( &_segment, _xyzfile );
       _segments.push_back(&_segment);
    }
    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
    _qmpackage->setLog( &_log );       
    _qmpackage->Initialize( &_package_options );
    // set the run dir 
    _qmpackage->setRunDir(".");
    Orbitals _orbitals;

   // do GWBSE   
   ExcitationEnergies( _qmpackage, _segments, &_orbitals);
   
 

   // only do this, if optimization is requested 
   if ( _do_optimize ){
   
      // optimization implies running DFT for displaced geometries
      _do_dft_input = true;
      _do_dft_run = true;
      _do_dft_parse = true;
      
      bool _converged = false;
      int _iteration = 0;
      // get total energy for this excited state
      double energy = GetTotalEnergy( &_orbitals, _spintype, _opt_state );
      cout << " energy of " << _spintype << " " << _opt_state << " is " << energy << " a.u. at iteration " << _iteration << endl;
      double energy_new;
      // now iterate
      while ( ! _converged ){ 
   
           vector <Segment* > _molecule;
           Segment _current_coordinates(0, "mol");
           // get current coordinates from Orbitals object
           Orbitals2Segment( &_current_coordinates , &_orbitals);
           _molecule.push_back( &_current_coordinates );
           
           // displace all atoms in each cartesian coordinate and get new energy
           // for each atom
           vector< Atom* > _atoms;
           vector< Atom* > ::iterator ait;
           _atoms = _current_coordinates.Atoms();
           
           if ( _iteration == 0 ){
              _natoms = _current_coordinates.Atoms().size();
              _force = ub::zero_matrix<double>(_natoms ,3);
              _force_old = ub::zero_matrix<double>(_natoms ,3);
              _xyz_shift = ub::zero_matrix<double>(_natoms ,3);
              _speed = ub::zero_matrix<double>(_natoms ,3);
              _current_xyz = ub::zero_matrix<double>(_natoms ,3);
              _old_xyz = ub::zero_matrix<double>(_natoms ,3);
           }
           
           
           int _i_atom = 0;
           
           for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
               // get its current coordinates
               vec _current_pos = (*ait)->getQMPos(); // in nm
               
               for ( int _i_cart = 0 ; _i_cart < 3; _i_cart++ ){
                   
                            // get displacement vector
                            vec _displaced(0, 0, 0);
                            if (_i_cart == 0) {
                                _displaced.setX(_displacement / 10); // x, _displacement in in Angstrom
                                _current_xyz(_i_atom,_i_cart ) = _current_pos.getX() * 18.897259886; // in Bohr
                            }
                            if (_i_cart == 1) {
                                _current_xyz(_i_atom,_i_cart ) = _current_pos.getY() * 18.897259886; // in Bohr
                                _displaced.setY(_displacement / 10); // y, _displacement in in Angstrom
                            }
                            if (_i_cart == 2) {
                                _current_xyz(_i_atom,_i_cart ) = _current_pos.getZ() * 18.897259886; // in Bohr
                                _displaced.setZ(_displacement / 10); // z, _displacement in in Angstrom
                            }
                            
                            // update the coordinate
                            vec _pos_displaced = _current_pos + _displaced;
                            (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                            // run DFT and GW-BSE for this geometry
                            ExcitationEnergies(_qmpackage, _molecule, &_orbitals);
                            
                            // get total energy for this excited state
                            double energy_displaced = GetTotalEnergy(&_orbitals, _spintype, _opt_state);
                            // calculate force and put into matrix
                            _force(_i_atom,_i_cart) = (energy - energy_displaced) / _displacement / 1.889725989; // force a.u./a.u.
                                                  
                            (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
               }
               
               cout << " Force on atom " << _i_atom << " " << setprecision(12) << _force(_i_atom, 0 ) << " : " << setprecision(12) << _force(_i_atom, 1 ) << " : " << _force(_i_atom, 2 ) << endl;
               _i_atom++;

           
        }
           
               ConjugateGradientStep( _iteration, _force, _force_old, _xyz_shift, _speed, _current_xyz, _old_xyz );
           
               cout << " Atoms are shifted by (in Ang): " << endl;
               
               for ( int _i_atoms = 0 ; _i_atoms < _natoms; _i_atoms++){
                   cout << "  atom " << _i_atoms << " : x " << _xyz_shift(_i_atoms,0)*0.529177249 << " y " << _xyz_shift(_i_atoms,1)*0.529177249  << " z " << _xyz_shift(_i_atoms,2)*0.529177249 << endl;
               }
           
                
           
           // put new coordinates into _segment
           _i_atom = 0;
           for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
               // put new coordinates (_current_xyz is in Bohr, segments in nm)
               vec _pos_displaced( _current_xyz(_i_atom,0)*0.052917725 , _current_xyz(_i_atom,1)*0.052917725, _current_xyz(_i_atom,2)*0.052917725 );
               (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment
               _i_atom++;
           }

           // for the updated geometry, get new reference energy
           ExcitationEnergies(_qmpackage, _molecule, &_orbitals);
           double energy_new  = GetTotalEnergy( &_orbitals, _spintype, _opt_state );
           cout << " energy of " << _spintype << " " << _opt_state << " is " << energy_new << " a.u. at iteration " << _iteration+1 << endl;
           
           double energy_delta = energy_new - energy;
           cout << "    energy change is "  << setprecision(12) << energy_delta << " a.u. or  " << energy_delta*27.21138386 << " eV " << endl;
           energy = energy_new;
           if ( std::abs(energy_delta*27.21138386) < _convergence) _converged = true;
           _iteration++;
       } // optimization while loop
       
   }// if optimization
  
           
            
       LOG(logDEBUG,_log) << "Saving data to " << _output_file << flush;
       std::ofstream ofs( ( _output_file).c_str() );
       boost::archive::binary_oarchive oa( ofs );


       
       oa << _orbitals;
       ofs.close();

     
     
     
     
    Property _summary; 
    Property *_job_output = &_summary.add("output","");
    votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
     
    //ofs (_output_file.c_str(), std::ofstream::out);
    //ofs << *_job_output;    
    //ofs.close();
    
    return true;
}




void Exciton::ConjugateGradientStep(int& _iteration, ub::matrix<double>& _force, ub::matrix<double>& _force_old, ub::matrix<double>& _xyz_shift, ub::matrix<double>& _speed, ub::matrix<double>& _current_xyz, ub::matrix<double>& _old_xyz){
    

    // stuff to prepare
    ub::vector<double> _total_force(3,0.0);
    double _force_norm = 0.0;
    for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
        for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
           _total_force( _i_cart ) += _force(_i_atom,_i_cart);
           _force_norm += _force(_i_atom,_i_cart)*_force(_i_atom,_i_cart);
        }
    }
    
    cout << " Total force " << _total_force( 0 ) << " " << _total_force( 1 ) << " " << _total_force( 2 ) << endl;
    cout << " Total force norm " << _force_norm;
            
            
    // remove CoM force
    for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
        for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
          _force(_i_atom,_i_cart) -= _total_force(_i_cart)/_natoms;
        }
    }
    
    
    cout << " Forces after CoM removal:" << endl;
    for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
        cout << " Force on  " << _i_atom << _force( _i_atom, 0 ) << " " << _force( _i_atom, 1 ) << " " << _force( _i_atom, 2 ) << endl;
    }
    
    // first trial point is just a plain step in the direction of the forces
    if ( _iteration == 0 ){
        _force_old = _force; // backup force
        _old_xyz = _current_xyz; // backup coordinates
        _speed = _force; // initialize 
        double _lambda0 = 0.2;
        double _factor = _lambda0/sqrt(_force_norm);
        _xyz_shift = _factor * _force;
        // new xyz
        _current_xyz += _xyz_shift;
    } else if  ( _iteration % 2 != 0 ) {
        // odd iterations (1,3,5,...) => minimize force along the line
        ub::matrix<double> _step = _xyz_shift; // backup copy of _xyz_shift
        double alpha1 = 0.0;
        double alpha2 = 0.0;
        for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
            for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
              alpha1 += _speed(_i_atom, _i_cart) * _force_old(_i_atom,_i_cart);
              alpha2 += _speed(_i_atom, _i_cart) * (_force(_i_atom,_i_cart) - _force_old(_i_atom,_i_cart));
            }  
        }
        double lambda = -alpha1/alpha2; // optimal lambda
        _xyz_shift = lambda * _step; // 
        // new coordinates
        _current_xyz = _old_xyz + _xyz_shift; 
        
    } else {
        // even iteration (2,4,6,...) => construct new direction and trial point
        ub::matrix<double> _step = _xyz_shift; // backup copy of _xyz_shift
        double alpha1 = 0.0;
        double alpha2 = 0.0;
        for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
            for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
              alpha1 += _speed(_i_atom, _i_cart) * _force_old(_i_atom,_i_cart);
              alpha2 += _speed(_i_atom, _i_cart) * (_force(_i_atom,_i_cart) - _force_old(_i_atom,_i_cart));
            }  
        }
        double lambda = -alpha1/alpha2; // optimal lambda
        _xyz_shift = lambda * _step; // 
        // new coordinates
        _current_xyz = _old_xyz + _xyz_shift; 
        
        // new step direction
        double _norm_old_force = 0.0;
        alpha1 = 0.0;
        for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
            for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
              _norm_old_force += _force_old(_i_atom,_i_cart) * _force_old(_i_atom,_i_cart);
              alpha1 += _force(_i_atom, _i_cart) * (_force(_i_atom,_i_cart) - _force_old(_i_atom,_i_cart));
            }  
        }
        double beta = _force_norm / _norm_old_force; // original conjugate-gradient
        // double beta = alpha1 / _norm_old_force; // Polak/Ribiere (?)
        _speed = beta*_speed + _force;
        
        // construct new step along new direction
        double _norm_h = 0.0;
        double _norm_old_step = 0.0;
        for ( int _i_atom = 0 ; _i_atom < _natoms; _i_atom++ ){
            for ( int _i_cart = 0; _i_cart < 3; _i_cart++ ){
              _norm_h += _speed(_i_atom,_i_cart) * _speed(_i_atom,_i_cart);
              _norm_old_step += _xyz_shift(_i_atom, _i_cart) * _xyz_shift(_i_atom,_i_cart);
            }  
        }
        double lambda0 = 2.0 * sqrt(_norm_old_step);
        beta = lambda0 / sqrt(_norm_h);
        
        // finally update stuff
        _old_xyz = _current_xyz;
        _current_xyz += beta*_speed;
        _xyz_shift = beta*_speed;
        _force_old = _force;
        
        //cout << " Higher iterations > 0 not implemented yet!" << endl;
        //exit(1);
        
    }

    
    
    
}


void Exciton::ExcitationEnergies(QMPackage* _qmpackage, vector<Segment*> _segments, Orbitals* _orbitals){


  if ( _do_dft_input ){
          _qmpackage->WriteInputFile( _segments );
       }
       if ( _do_dft_run ){
          _qmpackage->Run();
       }

      // parse DFT data, if required
      if ( _do_dft_parse ){
        LOG(logDEBUG,_log) << "Parsing DFT data " << _output_file << flush;
        _qmpackage->setOrbitalsFileName( _orbfile );
        int _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( _orbitals );
        _qmpackage->setLogFileName( _logfile );
        int _parse_log_status = _qmpackage->ParseLogFile( _orbitals );
     }

     // if no parsing of DFT data is requested, reload serialized orbitals object
     if ( !_do_dft_parse ){
         LOG(logDEBUG,_log) << "Loading serialized data from " << _output_file << flush;
         _orbitals->Load( _output_file );
     }

     if ( _do_gwbse ){
        _gwbse.setLogger(&_log);
        _gwbse.Initialize( &_gwbse_options );
        bool _evaluate = _gwbse.Evaluate( _orbitals );
        // std::cout << _log;
     }
     

    
    
}



double Exciton::GetTotalEnergy(Orbitals* _orbitals, string _spintype, int _opt_state){
    
    // total energy of the excited state
    double _total_energy ;
    double _omega;
    
    double _dft_energy = _orbitals->getQMEnergy();
    
    if ( _spintype == "singlet" ){
       _omega      = _orbitals->BSESingletEnergies()[_opt_state - 1];
    } else if ( _spintype == "triplet" ){
        _omega      = _orbitals->BSETripletEnergies()[_opt_state - 1];
    }

    
    // DFT total energy is stored in eV
    // singlet energies are stored in Ryd...
    return _total_energy = _dft_energy/27.21138386 + _omega*0.5; // a.u.
}


void Exciton::ReadXYZ(Segment* _segment, string filename){
 
    
               
                ifstream in;
                double x, y, z;
                int natoms, id;
                string label, type;
                vec pos;

                cout << " Reading molecular coordinates from " << _xyzfile << endl;
                in.open(_xyzfile.c_str(), ios::in);
                if (!in) throw runtime_error(string("Error reading coordinates from: ")
                        + _xyzfile);

                id = 1;
                while (in.good()) { // keep reading until end-of-file
                    in >> type;
                    in >> x;
                    in >> y;
                    in >> z;
                    if (in.eof()) break;
                    // cout << type << ":" << x << ":" << y << ":" << z << endl;

                    // creating atoms and adding them to the molecule
                    Atom *pAtom = new Atom(id++, type);
                    vec position(x / 10, y / 10, z / 10); // xyz has Angstrom, votca stores nm
                    pAtom->setPos(position);
                    pAtom->setQMPart(id, position);
                    pAtom->setElement(type);
                    _segment->AddAtom(pAtom);

                }
                in.close();
    
}

void Exciton::Orbitals2Segment(Segment* _segment, Orbitals* _orbitals){
    
            vector< QMAtom* > _atoms;
            vector< QMAtom* > ::iterator ait;
            _atoms = _orbitals->QMAtoms();
            
            string type;
            int id = 1;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                
                // Atom *pAtom = new Atom(id++, type);
                
                type = (*ait)->type;
                double x = (*ait)->x;
                double y = (*ait)->y;
                double z = (*ait)->z;
                Atom *pAtom = new Atom(id++, type);
                // cout << type << " " << x << " " << y << " " << z << endl;
                vec position(x/10, y/10, z/10); // xyz has Angstrom, votca stores nm
                pAtom->setPos(position);
                pAtom->setQMPart(id, position);
                pAtom->setElement(type);
                _segment->AddAtom(pAtom);
            }
            // cout << "??? " << _segment.Atoms().size() << endl;
}




}}


#endif
