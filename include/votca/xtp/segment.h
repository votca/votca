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


#ifndef __VOTCA_XTP_SEGMENT_H
#define	__VOTCA_XTP_SEGMENT_H

#include <votca/xtp/segmenttype.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/polarsite.h>
#include <votca/xtp/apolarsite.h>

class Topology;

namespace votca { namespace xtp {

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

    void             setHasState(bool yesno, int state);
    bool             hasState(int state);

    double           getOcc(int e_h_s_t);
    void             setOcc(double occ, int e_h_s_t);

    // state: -1 electron +1 hole +2 singlet +3 triplet
    
    void             setU_cC_nN(double dU, int state);
    void             setU_nC_nN(double dU, int state);
    void             setU_cN_cC(double dU, int state);
    void             setU_xX_nN(double dU, int state);
    void             setU_nX_nN(double dU, int state);
    void             setU_xN_xX(double dU, int state);
    const double    &getU_cC_nN(int state);
    const double    &getU_nC_nN(int state);
    const double    &getU_cN_cC(int state);
    const double    &getU_xX_nN(int state);
    const double    &getU_nX_nN(int state);
    const double    &getU_xN_xX(int state);
    double           getSiteEnergy(int state);

    double           getEMpoles(int state);
    void             setEMpoles(int state, double energy);
    bool             hasChrgState(int state) { return _hasChrgState[state+1]; }
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
    void             AddAPolarSite(APolarSite *pole);
    vector< Fragment* > &Fragments() { return _fragments; }
    vector < Atom* >    &Atoms() { return _atoms; }
    vector<PolarSite*>  &PolarSites() { return _polarSites; }
    vector<APolarSite*> &APolarSites() { return _apolarSites; }


    void Rigidify();

    void WritePDB(FILE *out, string tag1 = "Fragments", string tag2 = "MD");
    void WriteXYZ(FILE *out, bool useQMPos = true);

private:

    int         _id;
    string      _name;
    SegmentType *_typ;
    Topology    *_top;
    Molecule    *_mol;

    vector < Fragment* >    _fragments;
    vector < Atom* >        _atoms;
    vector < PolarSite* >   _polarSites;
    vector < APolarSite* >  _apolarSites;

    vec         _CoM;


    double _U_cC_nN_e;   // from ::EInternal     input     DEFAULT 0
    double _U_cC_nN_h;

    double _U_nC_nN_e;   // from ::EInternal     input     DEFAULT 0
    double _U_nC_nN_h;

    double _U_cN_cC_e;   // from ::EInternal     input     DEFAULT 0
    double _U_cN_cC_h;

    //double _ePolar_e;    // from ::EMultipole    output    DEFAULT 0
    //double _ePolar_h;

    double _occ_e;       // from ::KMC           output    DEFAULT 0
    double _occ_h;

    bool   _has_e;       // from ::EInternal     input     DEFAULT 0
    bool   _has_h;
    
    
    bool   _occ_s;      //state 3 = triplet
    bool   _occ_t;      // t:triplet s:singlet
    bool   _has_s;      // state 2 =. singlet
    
    
    bool   _has_t;      //Exciton Properties               DEFAULT 0
   
   
    
    double _U_xX_nN_s;
    double _U_xX_nN_t;
    
    double _U_nX_nN_s;   
    double _U_nX_nN_t;

    double _U_xN_xX_s;   
    double _U_xN_xX_t;
    
    //double _ePolar_s;
    //double _ePolar_t;


    vector< double > _eMpoles;
    //   +1(=> h)   e.static + pol. energy E(+1) - E(0)
    //   -1(=> e)   e.static + pol. energy E(-1) - E(0)
    //   +2(=> s)   e.static + pol. energy E(+2) - E(0)
    //   +3(=> t)   e.static + pol. energy E(+3) - E(0)
    vector<bool> _hasChrgState;

    map<int, vec> _intCoords;
    // qmid => position
    

};

}}

#endif	/* __VOTCA_XTP_SEGMENT_H */

