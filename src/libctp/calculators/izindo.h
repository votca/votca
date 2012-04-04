#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include <votca/ctp/qmpair2.h>
#include <votca/ctp/paircalculator2.h>
#include <votca/moo/mol_and_orb.h>

namespace votca { namespace ctp {

namespace MOO = votca::moo;

/**
* \brief Semi-empirical electronic coupling elements for QM pairs
*
* Semi-empirical (ZINDO) electronic coupling elements for all conjugated
* segments from the neighbout list. Requires molecular orbitals in GAUSSIAN
* format.
*
* Callname: izindo
*/

class IZindo : public ParallelPairCalculator
{
public:

    IZindo() {};
   ~IZindo() {};

    string  Identify() { return "IZindo"; }
    void    Initialize(Topology *top, Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    void    EvalPair(Topology *top, QMPair2 *pair, int slot);

    void    CTP2MOO2CTP(QMPair2 *pair, int slot, int state);
    void    CalculateJ(QMPair2 *pair);
    void    CleanUp();


private:

    MOO::mol_and_orb *_morb1;
    MOO::mol_and_orb *_morb2;
    MOO::basis_set   *_basis1;
    MOO::basis_set   *_basis2;
    MOO::orb         *_orb1;
    MOO::orb         *_orb2;
    MOO::fock        *_fock;

    vector<int> _torbNrs1;
    vector<int> _torbNrs2;

    // Fock::CalcJ(...) apparently uses
    // thread-unsafe containers and / or
    // (hidden) malloc()s. Restrict...
    Mutex       _FockPath;


    // Information on orbitals for segments
    map<string,string>          _seg_basisName;
    map<string,string>          _seg_orbFile;
    map<string,bool>            _seg_has_e;
    map<string,bool>            _seg_has_h;
    map< string, vector<int> >  _seg_torbs_e;
    map< string, vector<int> >  _seg_torbs_h;

};

void IZindo::CleanUp() {

    delete _morb1;
    delete _morb2;
    delete _basis1;
    delete _basis2;
    delete _orb1;
    delete _orb2;
    delete _fock;

    _torbNrs1.clear();
    _torbNrs2.clear();

}

void IZindo::Initialize(Topology *top, Property *options) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";

    /* ---- OPTIONS.XML Structure -----
     *
     * <izindo>
     *
     *      <orbitalsXML>ORBITALS.XML</orbitalsXML>
     *
     * </izindo>
     *
     */

    this->ParseOrbitalsXML(top, options);

}


void IZindo::ParseOrbitalsXML(Topology *top, Property *opt) {

    string key = "options.izindo";
    string orbitalsXML = opt->get(key+".orbitalsXML").as<string> ();
    cout << endl << "... ... Orbital data from " << orbitalsXML << ". ";

    Property alloc;
    load_property_from_xml(alloc, orbitalsXML.c_str());    

    /* --- ORBITALS.XML Structure ---
     *
     * <topology>
     *
     *     <molecules>
     *          <molecule>
     *          <name></name>
     *
     *          <segments>
     *
     *              <segment>
     *              <name></name>
     *
     *              <basisset></basisset>
     *              <orbitals></orbitals>
     *
     *              <torbital_e></torbital_e>
     *              <torbital_h></torbital_h>
     *              </segment>
     *
     *              <segment>
     *                  ...
     *
     */

    key = "topology.molecules.molecule";
    list<Property*> mols = alloc.Select(key);
    list<Property*> ::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); ++molit) {

        key = "segments.segment";
        list<Property*> segs = (*molit)->Select(key);
        list<Property*> ::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); ++segit) {

            string segName = (*segit)->get("name").as<string> ();
            string basisName = (*segit)->get("basisset").as<string> ();
            string orbFile = (*segit)->get("orbitals").as<string> ();

            bool has_e = false;
            bool has_h = false;
            vector<int> torbs_e;
            vector<int> torbs_h;

            if ( (*segit)->exists("torbital_e") ) {
                torbs_e = (*segit)->get("torbital_e").as< vector<int> > ();
                has_e = (torbs_e.size()) ? true : false;
            }
            if ( (*segit)->exists("torbital_h") ) {
                torbs_h = (*segit)->get("torbital_h").as< vector<int> > ();
                has_h = (torbs_h.size()) ? true : false;
            }

            _seg_basisName[segName] = basisName;
            _seg_orbFile[segName] = orbFile;
            _seg_has_e[segName] = has_e;
            _seg_has_h[segName] = has_h;
            _seg_torbs_e[segName] = torbs_e;
            _seg_torbs_h[segName] = torbs_h;

        }
    }
}


void IZindo::EvalPair(Topology *top, QMPair2 *qmpair, int slot) {

    this->LockCout();
    cout << "\r... ... Evaluating pair " << qmpair->getId()+1 << flush;
    this->UnlockCout();

    string segName1 = qmpair->Seg1()->getName();
    string segName2 = qmpair->Seg2()->getName();

    bool pair_has_e = false;
    bool pair_has_h = false;

    try {
        pair_has_e = _seg_has_e.at(segName1) && _seg_has_e.at(segName2);
        pair_has_h = _seg_has_h.at(segName1) && _seg_has_h.at(segName2);
    }
    catch (out_of_range) {
        this->LockCout();
        cout << endl << "... ... WARNING: No orbital information for pair ["
                     << segName1 << ", " << segName2 << "]. "
                     << "Skipping... " << endl;
        this->UnlockCout();

        return;
    }

    if (pair_has_e) {
        this->CTP2MOO2CTP(qmpair, slot, -1);
    }
    if (pair_has_h) {
        this->CTP2MOO2CTP(qmpair, slot, +1);
    }
}



void IZindo::CTP2MOO2CTP(QMPair2 *pair, int slot, int state) {

    // ++++++++++++++++++++++ //
    // Initialize MOO Objects //
    // ++++++++++++++++++++++ //

    Segment *seg1 = pair->Seg1PbCopy();
    Segment *seg2 = pair->Seg2PbCopy();

    SegmentType *type2 = seg2->getType();
    SegmentType *type1 = seg1->getType();

    string coordsFile1 = type1->getQMCoordsFile();
    string coordsFile2 = type2->getQMCoordsFile();

    string basisName1;
    string basisName2;
    string orbFile1;
    string orbFile2;
    vector<int> torbs1;
    vector<int> torbs2;
        
    basisName1 = _seg_basisName.at(type1->getName());
    basisName2 = _seg_basisName.at(type2->getName());

    orbFile1 = _seg_orbFile.at(type1->getName());
    orbFile2 = _seg_orbFile.at(type2->getName());

    torbs1 = (state == -1) ? _seg_torbs_e.at(type1->getName()) :
                             _seg_torbs_h.at(type1->getName());
    torbs2 = (state == -1) ? _seg_torbs_e.at(type2->getName()) :
                             _seg_torbs_h.at(type2->getName());
    
    //string basisName1 = type1->getBasisName();
    //string basisName2 = type2->getBasisName();
    //string orbFile1 = type1->getOrbitalsFile();
    //string orbFile2 = type2->getOrbitalsFile();
    //vector<int> torbs1 = type1->getTOrbNrs();
    //vector<int> torbs2 = type2->getTOrbNrs();
    
    MOO::mol_and_orb *morb1 = new MOO::mol_and_orb();
    MOO::mol_and_orb *morb2 = new MOO::mol_and_orb();
    MOO::basis_set *basis1 = new MOO::basis_set();
    MOO::basis_set *basis2 = new MOO::basis_set();
    MOO::orb *orb1   = new MOO::orb();
    MOO::orb *orb2   = new MOO::orb();
    MOO::fock *fock12   = new MOO::fock();

    vector<int> torbNrs1;
    vector<int> torbNrs2;

    // Define basis set //
    basis1->set_basis_set(basisName1);
    morb1->define_bs(basis1);
    basis2->set_basis_set(basisName2);
    morb2->define_bs(basis2);

    // Load QM coordinates //
    morb1->init(coordsFile1.c_str());
    morb2->init(coordsFile2.c_str());

    // Create orbitals //
    morb1->init_orbitals(*orb1, orbFile1.c_str());
    morb2->init_orbitals(*orb2, orbFile2.c_str());

    orb1->strip_orbitals(torbs1);
    orb2->strip_orbitals(torbs2);
    
    int frontier1 = torbs1.size();
    int frontier2 = torbs2.size();
    for (int i = 0; i < frontier1; i++) { torbNrs1.push_back(i); }
    for (int j = 0; j < frontier2; j++) { torbNrs2.push_back(j); }
    
    morb1->assign_orb(orb1);
    morb2->assign_orb(orb2);

    // Initialise Fock matrix //
    fock12->init(*morb1, *morb2);


    // ++++++++++++++++++++++++++++++++++++++++ //
    // Rotate + Translate to MD Frame: Mol&Orb1 //
    // ++++++++++++++++++++++++++++++++++++++++ //
    
    // Rotate + Translate QM frame to MD frame, fragment-wise
    vector< Fragment* > ::iterator fit;
    for (fit = seg1->Fragments().begin();
         fit < seg1->Fragments().end();
         fit++) {

        // Centers of maps, rotation matrix MD2QM
        Fragment *frag = *fit;
        vec CoMD = frag->getCoMD();
        vec CoQM = frag->getCoQM();
        matrix rotQM2MD = frag->getRotQM2MD();

        // Fill container with atom QM indices
        vector<int> atmIdcs;
        vector< Atom* > ::iterator ait;
        for (ait = frag->Atoms().begin();
             ait < frag->Atoms().end();
             ait++) {
             if ( (*ait)->HasQMPart() ) {
                 atmIdcs.push_back( (*ait)->getQMId()-1 );
             }
        }


        // Perform translation + rotation
        const double NM2Bohr= 10/0.529189379;
        morb1->rotate_someatoms(atmIdcs, rotQM2MD,
                                 CoMD*NM2Bohr, CoQM*NM2Bohr,
                                 morb1);

        // Rotate orbitals
        for (int i = 0; i < torbNrs1.size(); i++) {
            orb1->rotate_someatoms(atmIdcs, &rotQM2MD,
                                    morb1->getorb(i), i);
        }
    }
    
    //_morb1->write_pdb("morbs.pdb", "MOL", 1);

    
    // ++++++++++++++++++++++++++++++++++++++++ //
    // Rotate + Translate to MD Frame: Mol&Orb2 //
    // ++++++++++++++++++++++++++++++++++++++++ //

    for (fit = seg2->Fragments().begin();
            fit < seg2->Fragments().end();
            fit++) {

        // Centers of maps, rotation matrix MD2QM
        Fragment *frag = *fit;
        vec CoMD = frag->getCoMD();
        vec CoQM = frag->getCoQM();
        matrix rotQM2MD = frag->getRotQM2MD();

        // Fill container with atom QM indices
        vector<int> atmIdcs;
        vector< Atom* > ::iterator ait;
        for (ait = frag->Atoms().begin();
             ait < frag->Atoms().end();
             ait++) {
             if ( (*ait)->HasQMPart() ) {
                atmIdcs.push_back( (*ait)->getQMId()-1 );
             }
        }


        // Perform translation + rotation
        const double NM2Bohr= 10/0.529189379;
        morb2->rotate_someatoms(atmIdcs, rotQM2MD,
                                 CoMD*NM2Bohr, CoQM*NM2Bohr,
                                 morb2);

        // Rotate orbitals
        for (int i = 0; i < torbNrs2.size(); i++) {
            orb2->rotate_someatoms(atmIdcs, &rotQM2MD,
                                    morb2->getorb(i), i);
        }
    }

    // ++++++++++++++++++++++++++++ //
    // Calculate transfer integrals //
    // ++++++++++++++++++++++++++++ //

    
    vector<double> Js;
    std::pair< int, int > torb2torb;

    for (int i = 0; i < torbNrs1.size(); i++) {
    for (int j = 0; j < torbNrs2.size(); j++) {

        torb2torb.first  = torbNrs1[i];
        torb2torb.second = torbNrs2[j];

        _FockPath.Lock();
        double J = fock12->calcJ(torb2torb);
        _FockPath.Unlock();

        Js.push_back(J);
        
    }}

    pair->setJs(Js, state);
    pair->setIsPathCarrier(true, state);
    //cout << endl << "Jeff2 " << pair->getJeff2(state) << " ___ state " << state << endl;
    //cout << "path  " << pair->isPathCarrier(state) << endl;

    delete morb1;
    delete morb2;
    delete basis1;
    delete basis2;
    delete orb1;
    delete orb2;
    delete fock12;

    torbNrs1.clear();
    torbNrs2.clear();

}


void IZindo::CalculateJ(QMPair2 *pair) {

    Segment *seg1 = pair->Seg1PbCopy();
    Segment *seg2 = pair->Seg2PbCopy();

    // ++++++++++++++++++++++++++++++++++++++++ //
    // Rotate + Translate to MD Frame: Mol&Orb1 //
    // ++++++++++++++++++++++++++++++++++++++++ //

    // Rotate + Translate QM frame to MD frame, fragment-wise
    vector< Fragment* > ::iterator fit;
    for (fit = seg1->Fragments().begin();
         fit < seg1->Fragments().end();
         fit++) {

        // Centers of maps, rotation matrix MD2QM
        Fragment *frag = *fit;
        vec CoMD = frag->getCoMD();
        vec CoQM = frag->getCoQM();
        matrix rotQM2MD = frag->getRotQM2MD();

        // Fill container with atom QM indices
        vector<int> atmIdcs;
        vector< Atom* > ::iterator ait;
        for (ait = frag->Atoms().begin();
             ait < frag->Atoms().end();
             ait++) {
             if ( (*ait)->HasQMPart() ) {
                 atmIdcs.push_back( (*ait)->getQMId()-1 );
             }
        }


        // Perform translation + rotation
        const double NM2Bohr= 10/0.529189379;
        _morb1->rotate_someatoms(atmIdcs, rotQM2MD,
                                 CoMD*NM2Bohr, CoQM*NM2Bohr,
                                 _morb1);

        // Rotate orbitals
        for (int i = 0; i < this->_torbNrs1.size(); i++) {
            _orb1->rotate_someatoms(atmIdcs, &rotQM2MD,
                                    _morb1->getorb(i), i);
        }
    }

    //_morb1->write_pdb("morbs.pdb", "MOL", 1);


    // ++++++++++++++++++++++++++++++++++++++++ //
    // Rotate + Translate to MD Frame: Mol&Orb2 //
    // ++++++++++++++++++++++++++++++++++++++++ //

    for (fit = seg2->Fragments().begin();
            fit < seg2->Fragments().end();
            fit++) {

        // Centers of maps, rotation matrix MD2QM
        Fragment *frag = *fit;
        vec CoMD = frag->getCoMD();
        vec CoQM = frag->getCoQM();
        matrix rotQM2MD = frag->getRotQM2MD();

        // Fill container with atom QM indices
        vector<int> atmIdcs;
        vector< Atom* > ::iterator ait;
        for (ait = frag->Atoms().begin();
             ait < frag->Atoms().end();
             ait++) {
             if ( (*ait)->HasQMPart() ) {
                atmIdcs.push_back( (*ait)->getQMId()-1 );
             }
        }


        // Perform translation + rotation
        const double NM2Bohr= 10/0.529189379;
        _morb2->rotate_someatoms(atmIdcs, rotQM2MD,
                                 CoMD*NM2Bohr, CoQM*NM2Bohr,
                                 _morb2);

        // Rotate orbitals
        for (int i = 0; i < this->_torbNrs2.size(); i++) {
            _orb2->rotate_someatoms(atmIdcs, &rotQM2MD,
                                    _morb2->getorb(i), i);
        }
    }

    // ++++++++++++++++++++++++++++ //
    // Calculate transfer integrals //
    // ++++++++++++++++++++++++++++ //


    vector<double> Js;
    std::pair< int, int > torb2torb;

    for (int i = 0; i < _torbNrs1.size(); i++) {
    for (int j = 0; j < _torbNrs2.size(); j++) {

        torb2torb.first  = _torbNrs1[i];
        torb2torb.second = _torbNrs2[i];

        double J = _fock->calcJ(torb2torb);
        Js.push_back(J);
    }}

    assert(false); // pair->setJs(Js);
}

}}

#endif	/* _CALC_INTEGRALS_H */
