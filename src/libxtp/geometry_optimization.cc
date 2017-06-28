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

#include <votca/tools/linalg.h>
#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/forces.h>
#include <votca/xtp/bfgs-trm.h>

namespace votca {
    namespace xtp {

        void GeometryOptimization::Initialize(Property *options) {

            // checking if there is only one segment
            _nsegments = _segments.size();
            if (_nsegments > 1) throw runtime_error(string("\n Geometry optimization of more than 1 conjugated segment not supported. Stopping!"));
            
            _opt_state    = options->ifExistsReturnElseReturnDefault<int>(".state", 1);
            _spintype     = options->ifExistsReturnElseReturnDefault<string>(".spintype", "singlet");
            

            // pre-check optimizer method
            std::vector<string> choices={ "BFGS-TRM" };
            _optimizer         = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".optimizer.method", choices);
            _optimizer_options = options->get(".optimizer");

            // pre-check forces method
            choices={ "forward", "central" };
            _force_method  = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".forces.method", choices);
            // _displacement = options->get(".displacement").as<double> ();
            _force_options = options->get(".forces");
           
            // restart from saved history
            _restart_opt = options->ifExistsReturnElseReturnDefault<bool>(".restart", false);
            if (_restart_opt && !boost::filesystem::exists("optimization.restart")) {
                throw runtime_error(string("\n Restart requested but optimization.restart file not found!")); 
            };
            
            _natoms = _segments[0]->Atoms().size();
            _force = ub::zero_matrix<double>(_natoms, 3);
            _force_old = ub::zero_matrix<double>(_natoms, 3);
            _xyz_shift = ub::zero_matrix<double>(_natoms, 3);
            _current_xyz = ub::zero_matrix<double>(_natoms, 3);
            _old_xyz = ub::zero_matrix<double>(_natoms, 3);
            _trial_xyz = ub::zero_matrix<double>(_natoms, 3);
            _hessian = ub::zero_matrix<double>(3 * _natoms, 3 * _natoms);
            return;

        }

        void GeometryOptimization::Evaluate() {

            
            CTP_LOG(ctp::logINFO, *_pLog) << "Requested geometry optimization of excited state " << _spintype << " " << _opt_state << flush;
                        
            // get a force object
            Forces _force_engine(_gwbse_engine,_qmpackage, _segments, _orbitals);
            _force_engine.Initialize( &_force_options );
            _force_engine.setLog( _pLog );
            _force_engine.SetOptState( _opt_state );
            _force_engine.SetSpinType( _spintype  );

            
            // get the optimizer
            if (  _optimizer == "BFGS-TRM" ){
                
                BFGSTRM _bfgstrm(_gwbse_engine, _qmpackage, _segments, _orbitals, _force_engine);
                _bfgstrm.Initialize( &_optimizer_options );
                _bfgstrm.setLog( _pLog );
                _bfgstrm.Optimize( );
                
                
            }

            return;
        }

        void GeometryOptimization::Checkpoint(std::vector<ctp::Segment*>& _molecule) {

            // writing current coordinates and forces
            FILE *out = NULL;
            std::vector<ctp::Segment*> ::iterator sit;
            string FILENAME = "geometry_optimization.xyz";
            for (sit = _molecule.begin(); sit < _molecule.end(); ++sit) {
                if (_iteration == 0) {
                    out = fopen(FILENAME.c_str(), "w");
                } else {
                    out = fopen(FILENAME.c_str(), "a");
                }
                WriteIteration(out, (*sit));
            }

            fclose(out);
            return;
        }



        // write iteration

        void GeometryOptimization::WriteIteration(FILE *out, ctp::Segment* _segment) {

            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Atom* > _atoms = _segment->Atoms();
            int qmatoms = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                ++qmatoms;
            }

            fprintf(out, "%6d \n", qmatoms);
            fprintf(out, "at iteration %d \n", _iteration);

            int _i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                votca::tools::vec pos = (*ait)->getQMPos();
                string name = (*ait)->getElement();

                fprintf(out, "%2s %4.7f %4.7f %4.7f atom_vector %4.7f %4.7f %4.7f \n",
                        name.c_str(),
                        pos.getX()*10,
                        pos.getY()*10,
                        pos.getZ()*10,
                        _force(_i_atom, 0)*51.4220822067,
                        _force(_i_atom, 1)*51.4220822067,
                        _force(_i_atom, 2)*51.4220822067);
                _i_atom++;
                return;
            }

            return;
        }

    }
}
