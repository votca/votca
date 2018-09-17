/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
/// For earlier commit history see ctp commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_SEGMENT_H
#define	VOTCA_XTP_SEGMENT_H

#include <map>
#include <vector>

#include <votca/tools/vec.h>

namespace votca { namespace xtp {

class Atom;
class PolarSite;
class Fragment;
class SegmentType;
class Topology;
class Molecule;  

class Segment
{
public:

    Segment(int id, std::string name);
    Segment(Segment *stencil);
   ~Segment();

    const int       &getId() { return _id; }
    const std::string    &getName() { return _name; }

    const tools::vec       &getPos() const { return _CoM; }
    void             setPos(tools::vec pos) { _CoM = pos; }
    // This gets the center of mass from the MD positions of the atoms
    void             calcPos();
    void             TranslateBy(const tools::vec &shift);
    
    void            calcApproxSize();
    double          getApproxSize(){return _approxsize;}

    void             setHasState(bool yesno, int state);
    bool             hasState(int state);

    double           getOcc(int e_h_s_t);
    void             setOcc(double occ, int e_h_s_t);

    // state: -1 electron +1 hole +2 singlet +3 triplet
  
    /// Following notation can be observed in: 
    /// [1. Victor, R. et al. Microscopic Simulations of Charge Transport in Disordered Organic Semiconductors. J. Chem. Theory Comput. 7, 3335–3345 (2011).] 
    /// Labeling of the following methods follows the following semantics:
    /// U - Energy 
    /// n - neutral geometry
    /// N - neutral state
    /// c - charged geometry
    /// C - charged state
    /// x - excited geometry
    /// X - excited state

    /// UcC - UnN
    void             setU_cC_nN(double dU, int state);
    /// UnN - UnN
    void             setU_nC_nN(double dU, int state);
    /// UcN - UcC
    void             setU_cN_cC(double dU, int state);
    /// UxX - UnN
    void             setU_xX_nN(double dU, int state);
    /// UnX - UnN
    void             setU_nX_nN(double dU, int state);
    /// UxN - UxX
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
    void             setChrgStates(std::vector<bool> yesno) { _hasChrgState = yesno;}

    inline void      setTopology(Topology *container) { _top = container; }
    Topology        *getTopology() { return _top; }
    inline void      setMolecule(Molecule *container) { _mol = container; }
    Molecule        *getMolecule() { return _mol; }
    inline void      setType(SegmentType *type) { _typ = type; }
    SegmentType     *getType() { return _typ; }

    void             AddFragment( Fragment* fragment );
    void             AddAtom( Atom* atom );
    void             AddPolarSite(PolarSite *pole);
    std::vector< Fragment* > &Fragments() { return _fragments; }
    std::vector < Atom* >    &Atoms() { return _atoms; }
    std::vector<PolarSite*> &PolarSites() { return _polarSites; }


    void Rigidify();

private:

    int         _id;
    std::string      _name;
    SegmentType *_typ;
    Topology    *_top;
    Molecule    *_mol;

    std::vector < Fragment* >    _fragments;
    std::vector < Atom* >        _atoms;
    std::vector < PolarSite* >  _polarSites;

    tools::vec         _CoM;
    double   _approxsize;


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


    std::vector< double > _eMpoles;
    //   +1(=> h)   e.static + pol. energy E(+1) - E(0)
    //   -1(=> e)   e.static + pol. energy E(-1) - E(0)
    //   +2(=> s)   e.static + pol. energy E(+2) - E(0)
    //   +3(=> t)   e.static + pol. energy E(+3) - E(0)
    std::vector<bool> _hasChrgState;

    std::map<int, tools::vec> _intCoords;
    // qmid => position
    

};

}}

#endif	// VOTCA_XTP_SEGMENT_H 

