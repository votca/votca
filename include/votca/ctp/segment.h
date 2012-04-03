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

#include <votca/ctp/segmenttype.h>
#include <votca/ctp/fragment.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/polarsite.h>

class Topology;

namespace votca { namespace ctp {

class Molecule;  

 
class Segment
{
public:

    Segment(int id, string name);
    Segment(Segment *stencil);
   ~Segment();

    const int       &getId() { return _id; }
    const string    &getName() { return _name; }

    const vec       &getPos() const { return _CoM; }
    void             setPos(vec pos) { _CoM = pos; }
    void             calcPos();
    void             TranslateBy(const vec &shift);

    bool             hasState(int e_h);

    const double    &getOcc(int e_h);
    void             setOcc(int e_h, double occ);
    bool             hasOccProb() { return _hasOccProb; }

    const double    &getESite(int e_h);
    const double    &getESiteIntra(int e_h);
    void             setESiteIntra(int e_h, double energy);
    bool             hasEIntra() { return _hasESiteIntra; }

    const double    &getLambdaIntra(int state0, int state1);
    void             setLambdaIntra(int state0, int state1, double lambda);
    bool             hasLambda() { return _hasLambdas; }

    const double     getEMpoles(int e_h);
    void             setEMpoles(int e_h, double energy);
    bool             hasChrgState(int e_h) { return _hasChrgState[e_h+1]; }
    void             setChrgStates(vector<bool> yesno) { _hasChrgState = yesno;}

    inline void      setTopology(Topology *container) { _top = container; }
    Topology        *getTopology() { return _top; }
    inline void      setMolecule(Molecule *container) { _mol = container; }
    Molecule        *getMolecule() { return _mol; }
    inline void      setType(SegmentType *type) { _typ = type; }
    SegmentType     *getType() { return _typ; }

    void             AddFragment( Fragment* fragment );
    void             AddAtom( Atom* atom );
    void             AddPolarSite(PolarSite *pole);
    vector< Fragment* > &Fragments() { return _fragments; }
    vector < Atom* >    &Atoms() { return _atoms; }
    vector<PolarSite*>  &PolarSites() { return _polarSites; }


    void Rigidify();

    void WritePDB(FILE *out, string tag1 = "Fragments", string tag2 = "MD");

private:

    Topology    *_top;
    Molecule    *_mol;
    SegmentType *_typ;

    vector < Fragment* >    _fragments;
    vector < Atom* >        _atoms;
    vector < PolarSite* >   _polarSites;

    string      _name;
    int         _id;
    vec         _CoM;


    double _eSiteIntra_e;    // from ::Rates      input     DEFAULT 0
    double _eSiteIntra_h;

    double _lambdaIntra_e;   // from ::Rates      input     DEFAULT 0
    double _lambdaIntra_h;

    double _ePolar_e;        // from ::EMultipole output    DEFAULT 0
    double _ePolar_h;

    double _occ_e;           // from ::KMC        output    DEFAULT 0
    double _occ_h;

    bool   _has_e;
    bool   _has_h;


    map< int, double > _eSiteIntra;
    //   +1(=> h)     E(CAtion) - E(Neutral)
    //   -1(=> e)     E(Anion)  - E(Neutral)
    bool _hasESiteIntra;

    map< int, map < int, double > > _lambdasIntra;
    //   +1(=> h)    0   lambdaCN = UnC - UnN (dischrg)
    //   -1(=> e)    0   lambdaAN = UnA - UnN (dischrg)
    //    0         +1   lambdaNC = UcN - UcC (chrg)
    //    0         -1   lambdaNA = UaN - UaA (chrg)
    bool _hasLambdas;

    vector< double > _eMpoles;
    //   +1(=> h)   e.static + pol. energy E(+1) - E(0)
    //   -1(=> e)   e.static + pol. energy E(-1) - E(0)
    vector<bool> _hasChrgState;

    map< int,       double >      _occProb;
    //   +1(=> h)   occ.prob. for hole
    //   -1(=> e)   occ.prob. for electron
    bool _hasOccProb;


    map<int, vec> _intCoords;
    // qmid => position
    

};

}}

#endif	/* __VOTCA_CTP_SEGMENT_H */

