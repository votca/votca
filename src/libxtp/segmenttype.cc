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


#include <votca/xtp/segmenttype.h>

namespace votca { namespace xtp {

    /**

inline void get_orient(const vec & a, const vec & b, const vec & c, matrix & cg) {

    cg.set(0, 0, a.getX());
    cg.set(0, 1, b.getX());
    cg.set(0, 2, c.getX());
    cg.set(1, 0, a.getY());
    cg.set(1, 1, b.getY());
    cg.set(1, 2, c.getY());
    cg.set(2, 0, a.getZ());
    cg.set(2, 1, b.getZ());
    cg.set(2, 2, c.getZ());
    
}

SegmentType::~SegmentType() {
    
    cout << "calling the SegmentType destructor" << endl;

    vector < vector <int> >::iterator it_list_atoms_monomer = _list_atoms_monomer.begin();

    cout << "start clearing the list of atom monomers" << endl;
    for (; it_list_atoms_monomer != _list_atoms_monomer.end(); ++it_list_atoms_monomer) {
        it_list_atoms_monomer->clear();
    }
    
    cout << "start clearing the list of  monomers" << endl;
    _list_atoms_monomer.clear();
    
    cout << "start clearing the list of com monomers" << endl;
    _list_coms_monomer.clear();
    
    cout << "start clearing the list of ors monomers" <<endl;
    _list_ors_monomer.clear();
 
}

// Constructor
SegmentType::SegmentType(const char * namecoord, const char * nameorb,
        const char * nameneutr, const char * namecrg, string & basisset,
        const vector < int>& transorbs, const unsigned int &id,
        string name, vector < vector < int > > list_atoms_monomer,
        vector < vector < double > > list_weights_monomer) {

    _id = id;
    _name = name;
    vector < vector < int > > ::iterator it_mon;
    vector < vector < double > > ::iterator it_mon_weights;

    if (list_atoms_monomer.size() != list_weights_monomer.size())
        throw invalid_argument("number of atoms and weights do not match for charge unit <" + name + ">");

    int count = 0;
    for (it_mon = list_atoms_monomer.begin(), it_mon_weights = list_weights_monomer.begin();
            it_mon < list_atoms_monomer.end(); ++it_mon, ++it_mon_weights) {
        
        vector <int> list_mon;
        _list_atoms_monomer.push_back(list_mon);
        
        vector <int> ::iterator it_at;
        vector <double> ::iterator it_weight;
        vec com(0., 0., 0.);
        matrix m(0.);
        matrix orient;
        vec xprime, yprime, zprime;
        if (it_mon->size() != it_mon_weights->size())
            throw invalid_argument("number of atoms and weights do not match for charge unit <" + name + ">");
        
        double sum_weights = 0;
        for (it_at = it_mon->begin(), it_weight = it_mon_weights->begin();
                it_at != it_mon->end(); ++it_at, ++it_weight) {
            com = com + _intcoords.GetPos(*it_at)*(*it_weight);
            (_list_atoms_monomer[count]).push_back(*it_at);
            sum_weights += *it_weight;
        }

        com /= sum_weights;
        _list_coms_monomer.push_back(com);

        if (it_mon-> size() > 3) {
            for (it_at = it_mon->begin(), it_weight = it_mon_weights->begin();
                    it_at != it_mon->end(); ++it_at, ++it_weight) {
                if (*it_weight == 0) continue;
                vec v = _intcoords.GetPos(*it_at) - com;

                cout << "Adding atom " << *it_at << " to monomer " << count << endl;
                m[0][0] +=  v.getX() * v.getX();
                m[0][1] +=  v.getX() * v.getY();
                m[0][2] +=  v.getX() * v.getZ();
                m[1][1] +=  v.getY() * v.getY();
                m[1][2] +=  v.getY() * v.getZ();
                m[2][2] +=  v.getZ() * v.getZ();
            }

            m[1][0] = m[0][1];
            m[2][0] = m[0][2];
            m[2][1] = m[1][2];

            matrix::eigensystem_t es;
            m.SolveEigensystem(es);
            zprime = es.eigenvecs[0];


            xprime = _intcoords.GetPos((list_atoms_monomer[count])[1]) - _intcoords.GetPos((list_atoms_monomer[count])[0]);
            xprime.normalize(); //these 2 lines useless
            zprime.normalize();

            cout << "xaxis : " << xprime << endl;

            vec w = _intcoords.GetPos((list_atoms_monomer[count])[2]) - _intcoords.GetPos((list_atoms_monomer[count])[0]);
            w.normalize();

            cout << "CN bond: " << w << endl;

            if ((xprime^w) * zprime < 0) {
                zprime = vec(0., 0., 0.) - zprime;

                cout << "flipping the zaxis" << endl;
            }

            vec yprime = (zprime) ^ (xprime);
            xprime = yprime ^ (zprime);
            xprime.normalize();
            yprime.normalize();
            zprime.normalize();

            get_orient(xprime, yprime, zprime, orient);
            orient.Transpose();
        } else {
            orient.UnitMatrix();
        }
        _list_ors_monomer.push_back(orient);

        cout << "Making segment of type: " << _id << " initializing monomer " << count << " with normal "
                << zprime << " and plane: " << xprime << endl;
        count++;
    }

}

/// this will move each rigid fragment into positions[i] by rotating 
/// corresponding to two given vectors 
/// we assume that the pointer is to a suitable molecule... 

void SegmentType::rotate_fragments(vector < vec >::iterator it_pos, vector < vec >::iterator it_norm,
        vector <vec >::iterator it_plan, mol_and_orb * rotated_molecule) {
    vector < vector <int> >::iterator it_monomers;
    vector < vec >::iterator it_coms = _list_coms_monomer.begin();
    vector <matrix>::iterator it_ors = _list_ors_monomer.begin();

    int count = 0;
    vec old_norm;
    for (it_monomers = _list_atoms_monomer.begin(); it_monomers != _list_atoms_monomer.end();
            ++it_monomers, ++it_pos, ++it_norm, ++it_plan, ++it_coms, ++it_ors) {

        
        // we need to define the rotation matrix for the beads
        vec zprime = (*it_norm);
        vec yprime = zprime ^ (*it_plan);
        vec xprime = yprime ^ (zprime);

        xprime = xprime.normalize();
        yprime = yprime.normalize();

        vec zprime_orb = zprime;
        if (count > 0) {
            // check if the orbitals need "flipping"
            if (old_norm * (zprime) > 0.) zprime_orb = vec(0., 0., 0.) - zprime;
        }

        matrix orient_pos;
        matrix orient_orb;

        get_orient(xprime, yprime, zprime, orient_pos);
        get_orient(xprime, yprime, zprime_orb, orient_orb);

        matrix Orient_Pos = orient_pos * (*it_ors);
        matrix Orient_Orb = orient_orb * (*it_ors);

        cout << "NUMBER OF atoms in the the fragment: " << it_monomers->size() << "displacement" << *it_pos <<
                " internal center of nassL " << *it_coms << endl;
        cout << " type of SegmentType: " << _id << endl;
        cout << "normal vector: " << zprime << " and for orbitals: " << zprime_orb << " and plane vector: " << xprime << endl;

        _intcoords.rotate_someatoms(*it_monomers, Orient_Pos, *it_pos, *it_coms, rotated_molecule);

        old_norm = zprime;
        count++;
    }
}
**/

}}
