#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include <votca/ctp/qmpair2.h>
#include <votca/ctp/paircalculator2.h>
#include <votca/moo/mol_and_orb.h>

namespace votca { namespace ctp {

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
    void    EvalPair(Topology *top, QMPair2 *pair, int slot);

    void    CTP2MOO2CTP(QMPair2 *pair, int slot);
    void    CalculateJ(QMPair2 *pair);
    void    CleanUp();


private:

    mol_and_orb *_morb1;
    mol_and_orb *_morb2;
    basis_set   *_basis1;
    basis_set   *_basis2;
    orb         *_orb1;
    orb         *_orb2;
    fock        *_fock;

    vector<int> _torbNrs1;
    vector<int> _torbNrs2;

    // Fock::CalcJ(...) apparently uses
    // thread-unsafe containers and / or
    // (hidden) malloc()s. Restrict...
    Mutex _FockPath;

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
}


void IZindo::EvalPair(Topology *top, QMPair2 *qmpair, int slot) {

    this->LockCout();
    cout << "\r... ... Evaluating pair " << qmpair->getId()+1 << flush;
    this->UnlockCout();


    this->CTP2MOO2CTP(qmpair, slot);

}



void IZindo::CTP2MOO2CTP(QMPair2 *pair, int slot) {

    // ++++++++++++++++++++++ //
    // Initialize MOO Objects //
    // ++++++++++++++++++++++ //

    Segment *seg1 = pair->Seg1PbCopy();
    Segment *seg2 = pair->Seg2PbCopy();

    SegmentType *type2 = seg2->getType();
    SegmentType *type1 = seg1->getType();

    string coordsFile1 = type1->getQMCoordsFile();
    string coordsFile2 = type2->getQMCoordsFile();
    string basisName1 = type1->getBasisName();
    string basisName2 = type2->getBasisName();
    string orbFile1 = type1->getOrbitalsFile();
    string orbFile2 = type2->getOrbitalsFile();
    vector<int> torbs1 = type1->getTOrbNrs();
    vector<int> torbs2 = type2->getTOrbNrs();

    
    mol_and_orb *morb1 = new mol_and_orb();
    mol_and_orb *morb2 = new mol_and_orb();
    basis_set *basis1 = new basis_set();
    basis_set *basis2 = new basis_set();
    orb *orb1   = new orb();
    orb *orb2   = new orb();
    fock *fock12   = new fock();

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
    
    morb1->write_pdb("morbs.pdb", "MOL", 0);

    
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
    
    morb2->write_pdb("morbs.pdb", "MOL", 1);

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

    pair->setJs(Js);

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

    pair->setJs(Js);
}

}}

#endif	/* _CALC_INTEGRALS_H */
