/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_CTP_SEGMENT_H
#define	__VOTCA_CTP_SEGMENT_H

#include<votca/ctp/fragment.h>
#include<votca/ctp/atom.h>

class Topology;

namespace votca { namespace ctp {


class Molecule;  

 
    
/**
    \brief Conjugated segment. One conjugated segment contains several rigid fragments.

 * Apart from the position it has a vector with pointers to all fragments 
 * which belong to it. Substitutes QMBead of ctp and crgunittype of moo
 */
class Segment
{
public:
    /// Default constructor
    Segment(int id, string name);
    /// Default destructor
   ~Segment();

    const int       &getId();
    const string    &getName();
    const double    &getOcc() { return _occ; }
    void             setOcc(double occ) { _occ = occ; }



    void AddFragment( Fragment* fragment );
    void AddAtom( Atom* atom );

    vector< Fragment* > &Fragments() { return _fragments; }
    vector < Atom* > &Atoms() { return _atoms; }

    inline void setTopology(Topology *container) { _top = container; }
    inline void setMolecule(Molecule *container) { _mol = container; }

    Topology *getTopology() { return _top; }
    Molecule *getMolecule() { return _mol; }

    void         calcPos();
    void         setPos(vec pos) { _CoM = pos; }
    const vec   &getPos() const { return _CoM; }

    void WritePDB(FILE *out);

private:

    Topology    *_top;
    Molecule    *_mol;

    vector < Fragment* >    _fragments;
    vector < Atom* >        _atoms;

    string      _name;
    int         _id;
    vec         _CoM;
    double      _occ;
    

};

}}

#endif	/* __VOTCA_CTP_SEGMENT_H */

