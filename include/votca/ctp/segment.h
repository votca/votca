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

namespace votca { namespace ctp {

class Molecule;   
class Atom;   
    
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
   /**
     * get the id of the segment
     * \return segment id
     */
    const int &getId();
     /**
     * get the name of the segment
     * \return segment name
     */
    const string &getName();


    /// Adds a pointer to a fragment belonging to this segment
    void AddFragment( Fragment* fragment );
    vector< Fragment* > &Fragments() { return _fragments; }

    /// Adds a pointer to an atom belonging to this segment
    void AddAtom( Atom* atom );
    int NumberOfAtoms() { return _atoms.size(); }

private:
    string _name;
    /// Conjugated segment ID
    int _id;  
    /// position of a segment
    vec _pos;
    /// List of rigid fragments which belong to this segment
    vector < Fragment* > _fragments;
    /// List of atoms which belong to this segment
    vector < Atom* > _atoms;
   /// Molecule this Segment belongs to
    Molecule *_molecule;
    /// Name of the conjugated segment       
};

}}

#endif	/* __VOTCA_CTP_SEGMENT_H */

