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

class IZindo : public PairCalculator2
{
public:

    IZindo() {};
   ~IZindo() {};

    string  Identify() { return "IZindo"; }
    void    EvaluatePair(Topology *top, QMPair2 *pair);
    
    void    TranslateCTP2MOO(QMPair2 *pair);
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

void IZindo::EvaluatePair(Topology *top, QMPair2 *pair) {

    cout << "\r... ... Evaluating pair " << pair->getId()+1 << flush;

    this->TranslateCTP2MOO(pair);
    this->CalculateJ(pair);
    this->CleanUp();

}




void IZindo::TranslateCTP2MOO(QMPair2 *pair) {

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

    _morb1 = new mol_and_orb();
    _morb2 = new mol_and_orb();
    _basis1 = new basis_set();
    _basis2 = new basis_set();
    _orb1   = new orb();
    _orb2   = new orb();
    _fock   = new fock();


    // Define basis set //
    _basis1->set_basis_set(basisName1);
    _morb1->define_bs(_basis1);
    _basis2->set_basis_set(basisName2);
    _morb2->define_bs(_basis2);


    // Load QM coordinates //
    _morb1->init(coordsFile1.c_str());
    _morb2->init(coordsFile2.c_str());    

    // Create orbitals //
    _morb1->init_orbitals(*_orb1, orbFile1.c_str());
    _morb2->init_orbitals(*_orb2, orbFile2.c_str());

    _orb1->strip_orbitals(torbs1);
    _orb2->strip_orbitals(torbs2);
    
    int frontier1 = torbs1.size();
    int frontier2 = torbs2.size();
    for (int i = 0; i < frontier1; i++) { _torbNrs1.push_back(i); }
    for (int j = 0; j < frontier2; j++) { _torbNrs2.push_back(j); }

    _morb1->assign_orb(_orb1);
    _morb2->assign_orb(_orb2);   


    // Initialise Fock matrix //
    _fock->init(*_morb1, *_morb2);

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

