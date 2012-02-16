#ifndef _CALC_INTEGRALS_H
#define	_CALC_INTEGRALS_H

#include <votca/ctp/parallelpaircalc.h>
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
    void    PrepareFrame(Topology *top);
    void    EvalPair(Topology *top, QMPair2 *qmpair, int slot);

    void    TranslateCTP2MOO(QMPair2 *pair, int slot);
    void    CalculateJ(QMPair2 *pair, int slot);
    void    CleanUp(int slot);


private:

    vector <mol_and_orb *> _morb1;
    vector <mol_and_orb *> _morb2;
    vector <basis_set *>   _basis1;
    vector <basis_set *>   _basis2;
    vector <orb *>         _orb1;
    vector <orb *>         _orb2;
    vector <fock *>        _fock;

    vector <vector<int> >  _torbNrs1;
    vector <vector<int> >  _torbNrs2;

};

void IZindo::CleanUp(int slot) {

    delete _morb1[slot];
    delete _morb2[slot];
    delete _basis1[slot];
    delete _basis2[slot];
    delete _orb1[slot];
    delete _orb2[slot];
    delete _fock[slot];

    _torbNrs1[slot].clear();
    _torbNrs2[slot].clear();

}

void IZindo::Initialize(Topology *top, Property *options) {

    _nThreads = 4;
    cout << endl << "... ... Initialize with " << _nThreads << " threads ";
}


void IZindo::PrepareFrame(Topology *top) {

    for (int i = 0; i < _nThreads; i++) {

        _morb1.push_back(NULL);
        _morb2.push_back(NULL);
        _basis1.push_back(NULL);
        _basis2.push_back(NULL);
        _orb1.push_back(NULL);
        _orb2.push_back(NULL);
        _fock.push_back(NULL);
        
        vector<int> torbNrs1;
        vector<int> torbNrs2;

        _torbNrs1.push_back(torbNrs1);
        _torbNrs2.push_back(torbNrs2);
    }
}


void IZindo::EvalPair(Topology *top, QMPair2 *qmpair, int slot) {

    LockCout();
    cout << "\r... ... Evaluating pair " << qmpair->getId()+1 << flush;
    UnlockCout();

    TranslateCTP2MOO(qmpair, slot);
    CalculateJ(qmpair, slot);
    CleanUp(slot);

}




void IZindo::TranslateCTP2MOO(QMPair2 *pair, int slot) {

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

    _morb1[slot]  = new mol_and_orb();
    _morb2[slot]  = new mol_and_orb();
    _basis1[slot] = new basis_set();
    _basis2[slot] = new basis_set();
    _orb1[slot]   = new orb();
    _orb2[slot]   = new orb();
    _fock[slot]   = new fock();


    // Define basis set //
    _basis1[slot]->set_basis_set(basisName1);
    _morb1[slot]->define_bs(_basis1[slot]);
    _basis2[slot]->set_basis_set(basisName2);
    _morb2[slot]->define_bs(_basis2[slot]);


    // Load QM coordinates //
    _morb1[slot]->init(coordsFile1.c_str());
    _morb2[slot]->init(coordsFile2.c_str());

    // Create orbitals //
    _morb1[slot]->init_orbitals(*_orb1[slot], orbFile1.c_str());
    _morb2[slot]->init_orbitals(*_orb2[slot], orbFile2.c_str());

    _orb1[slot]->strip_orbitals(torbs1);
    _orb2[slot]->strip_orbitals(torbs2);
    
    int frontier1 = torbs1.size();
    int frontier2 = torbs2.size();
    for (int i = 0; i < frontier1; i++) { _torbNrs1[slot].push_back(i); }
    for (int j = 0; j < frontier2; j++) { _torbNrs2[slot].push_back(j); }

    _morb1[slot]->assign_orb(_orb1[slot]);
    _morb2[slot]->assign_orb(_orb2[slot]);


    // Initialise Fock matrix //
    _fock[slot]->init(*_morb1[slot], *_morb2[slot]);

}


void IZindo::CalculateJ(QMPair2 *pair, int slot) {

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
        _morb1[slot]->rotate_someatoms(atmIdcs, rotQM2MD,
                                 CoMD*NM2Bohr, CoQM*NM2Bohr,
                                 _morb1[slot]);

        // Rotate orbitals
        for (int i = 0; i < this->_torbNrs1[slot].size(); i++) {
            _orb1[slot]->rotate_someatoms(atmIdcs, &rotQM2MD,
                                    _morb1[slot]->getorb(i), i);
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
        _morb2[slot]->rotate_someatoms(atmIdcs, rotQM2MD,
                                 CoMD*NM2Bohr, CoQM*NM2Bohr,
                                 _morb2[slot]);

        // Rotate orbitals
        for (int i = 0; i < this->_torbNrs2[slot].size(); i++) {
            _orb2[slot]->rotate_someatoms(atmIdcs, &rotQM2MD,
                                    _morb2[slot]->getorb(i), i);
        }        
    }

    // ++++++++++++++++++++++++++++ //
    // Calculate transfer integrals //
    // ++++++++++++++++++++++++++++ //


    vector<double> Js;
    std::pair< int, int > torb2torb;

    for (int i = 0; i < _torbNrs1[slot].size(); i++) {
    for (int j = 0; j < _torbNrs2[slot].size(); j++) {

        torb2torb.first  = _torbNrs1[slot][i];
        torb2torb.second = _torbNrs2[slot][i];

        double J = _fock[slot]->calcJ(torb2torb);
        Js.push_back(J);
    }}

    pair->setJs(Js);
}

}}

#endif	/* _CALC_INTEGRALS_H */

