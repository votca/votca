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
#include <votca/xtp/forces.h>

namespace votca {
    namespace xtp {

        void Forces::Initialize(Property *options) {

            // checking if there is only one segment
            _nsegments = _segments.size();
            if (_nsegments > 1) throw runtime_error(string("\n Force calculation for more than 1 conjugated segment not supported. Stopping!"));

            // pre-check forces method
            std::vector<string> choices = {"forward", "central"};
            _force_method = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(".method", choices);

            if ((_force_method == "forward") || (_force_method == "central")) {
                _displacement = options->ifExistsReturnElseReturnDefault<double>(".displacement", 0.001); // Angstrom
            }

            _natoms = _segments[0]->Atoms().size();
            _forces = ub::zero_matrix<double>(_natoms, 3);

            return;

        }

        void Forces::Calculate( const double& energy ) {


            cout << " Numerical forces " << _force_method << " with displacement " << _displacement << endl;

            //backup current coordinates (WHY?)
            std::vector <ctp::Segment* > _molecule;
            ctp::Segment _current_coordinates(0, "mol");
            Orbitals2Segment(&_current_coordinates, _orbitals);
            _molecule.push_back(&_current_coordinates);

            // displace all atoms in each Cartesian coordinate and get new energy
            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            _atoms = _current_coordinates.Atoms();

            int _i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                // Calculate Force on this atom
                ub::matrix_range< ub::matrix<double> > _force = ub::subrange(_forces, _i_atom, _i_atom + 1, 0, 3);
                if (_force_method == "forward") NumForceForward(energy, ait, _force, _molecule);

                _i_atom++;
            }
            return;
        }

        /* Calculate forces on an atom numerically by forward differences */
        void Forces::NumForceForward(double energy, std::vector< ctp::Atom* > ::iterator ait, ub::matrix_range< ub::matrix<double> >& _force,
                std::vector<ctp::Segment*> _molecule) {

            // get this atoms's current coordinates
            vec _current_pos = (*ait)->getQMPos(); // in nm
            
            for (int _i_cart = 0; _i_cart < 3; _i_cart++) {

                // get displacement std::vector
                vec _displaced(0, 0, 0);
                if (_i_cart == 0) {
                    _displaced.setX(_displacement/10.0 ); // x, _displacement in Angstrom, now in nm
                }
                if (_i_cart == 1) {
                    _displaced.setY(_displacement/10.0 ); // y, _displacement in in Angstrom, now in nm
                }
                if (_i_cart == 2) {
                    _displaced.setZ(_displacement/10.0 ); // z, _displacement in in Angstrom, now in nm
                }

                // update the coordinate
                vec _pos_displaced = _current_pos + _displaced;
                
                (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                // run DFT and GW-BSE for this geometry
                _gwbse_engine.ExcitationEnergies(_qmpackage, _molecule, _orbitals);

                // get total energy for this excited state
                double energy_displaced = _orbitals->GetTotalEnergy(_spin_type, _opt_state);
                
                // calculate force and put into matrix
                _force(0,_i_cart) = (energy - energy_displaced) / (_displacement*votca::tools::conv::ang2bohr); // force a.u./a.u.
                (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
            } // Cartesian directions
            return;
        }
        



void Forces::Orbitals2Segment(ctp::Segment* _segment, Orbitals* _orbitals) {

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
        vec position(x / 10, y / 10, z / 10); // xyz has Angstrom, votca stores nm
        pAtom->setPos(position);
        pAtom->setQMPart(id, position);
        pAtom->setElement(type);
        _segment->AddAtom(pAtom);
    }
    return;
}

}
}
