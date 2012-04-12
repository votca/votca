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

#ifndef __VOTCA_CTP_TOPOLOGY_H
#define	__VOTCA_CTP_TOPOLOGY_H

#include <votca/tools/property.h>

#include <votca/csg/boundarycondition.h>
#include <votca/csg/openbox.h>
#include <votca/csg/orthorhombicbox.h>
#include <votca/csg/triclinicbox.h>

#include <votca/ctp/polarsite.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/fragment.h>
#include <votca/ctp/segmenttype.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/molecule.h>

#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmnblist.h>

#include <votca/moo/jcalc.h>
#include <votca/moo/mol_and_orb.h>

namespace CSG = votca::csg;


namespace votca { namespace ctp {

/**
 * \brief Container for molecules, conjugated segments, rigid fragments,
 * and atoms.
*/
class Topology 
{
public:

    Topology();
   ~Topology();


    // Population: Molecules, segments, fragments, atoms

    Molecule    *AddMolecule (string molecule_name);
    Segment     *AddSegment  (string segment_name);
    Atom        *AddAtom     (string atom_name);
    Fragment    *AddFragment (string fragment_name);
    PolarSite   *AddPolarSite(string siteName);
    SegmentType *AddSegmentType (string typeName);

    Molecule    *getMolecule(int id) { return _molecules[id-1]; }
    Segment     *getSegment(int id)  { return _segments[id-1]; }
    Fragment    *getFragment(int id) { return _fragments[id-1]; }
    Atom        *getAtom(int id)     { return _atoms[id-1]; }
    PolarSite   *getPolarSite(int id) { return _polarSites[id-1]; }
    SegmentType *getSegmentType(int id) { return _segmentTypes[id-1]; }

    vector< Atom* >         &Atoms() { return _atoms; }
    vector< Fragment* >     &Fragments() { return _fragments; }
    vector< Segment* >      &Segments() { return _segments; }
    vector< Molecule* >     &Molecules() { return _molecules; }
    vector< PolarSite* >    &PolarSites() { return _polarSites; }
    vector< SegmentType* >  &SegmentTypes() { return _segmentTypes; }

    bool        Rigidify();
    void        setCanRigidify(bool yesno) { _canRigidify = yesno; }
    const bool  canRigidify() { return _canRigidify; }
    const bool  isRigid() { return _isRigid; }
    void        setIsEStatified(bool yesno) { _isEStatified = yesno; }
    const bool  isEStatified() { return _isEStatified; }



    // Periodic boundary: Can be 'open', 'orthorhombic', 'triclinic'

    vec              PbShortestConnect(const vec &r1, const vec &r2) const;
    const matrix    &getBox() { return _bc->getBox(); }
    double           BoxVolume() { _bc->BoxVolume(); }
    void             setBox(const matrix &box,
                            CSG::BoundaryCondition::eBoxtype boxtype =
                            CSG::BoundaryCondition::typeAuto);

    QMNBList       &NBList() { return _nblist; }

    // Trajectory meta data: step number, time, frame (= Db ID)

    const int        getStep() { return _step; }
    void             setStep(int step) { _step = step; }
    const double     getTime() { return _time; }
    void             setTime(double time) { _time = time; }

    int              getDatabaseId() { return _db_id; };
    void             setDatabaseId(int id) { _db_id = id; }
    void             CleanUp();

    void             PrintInfo(ostream &out);
    void             PrintInfo(FILE *out);
    void             WritePDB(FILE *out, string tag = "segments");

   
protected:

    vector < Molecule* >    _molecules;
    vector < Segment* >     _segments;
    vector < Fragment* >    _fragments;
    vector < Atom* >        _atoms;
    vector < PolarSite* >   _polarSites;
    vector < SegmentType* > _segmentTypes;

    QMNBList               _nblist;

    CSG::BoundaryCondition *_bc;
    bool                    _hasPb;

    bool                    _canRigidify;
    bool                    _isRigid;
    bool                    _isEStatified;

    double _time;
    int    _step;
    int    _db_id;


    CSG::BoundaryCondition::eBoxtype
    AutoDetectBoxType(const matrix &box);

};

}}

#endif	/* __VOTCA_CTP_TOPOLOGY_H */

