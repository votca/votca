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

    string Identify() { return "IZindo"; }    
    void EvaluatePair(Topology *top, QMPair2 *pair);
    
    
    void Translate(Segment *seg1, Segment *seg2);
    void CleanUp();





private:

    mol_and_orb *_morb1;
    mol_and_orb *_morb2;
    basis_set   *_basis1;
    basis_set   *_basis2;
    orb         *_orb1;
    orb         *_orb2;

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

    _torbNrs1.clear();
    _torbNrs2.clear();

}

void IZindo::EvaluatePair(Topology *top, QMPair2 *pair) {

    cout << "\r... ... Evaluating pair " << pair->getId()+1 << flush;

    // TODO Ensure PB ghost initialised correctly
    Segment *seg1 = pair->Seg1PbCopy();
    Segment *seg2 = pair->Seg2PbCopy();

    seg1->Rigidify();
    seg2->Rigidify();

    this->Translate(seg1, seg2);



    this->CleanUp();

//    CrgUnit *crg1 = pair->Crg1PBCCopy();
//    CrgUnit *crg2 = pair->Crg2PBCCopy();
//    vector <double> Js = top->GetJCalc().CalcJ(*crg1, *crg2);
//    pair->setJs(Js);


}




void IZindo::Translate(Segment *seg1, Segment *seg2) {

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

    _basis1->set_basis_set(basisName1);
    _morb1->define_bs(_basis1);
    _basis2->set_basis_set(basisName2);
    _morb2->define_bs(_basis2);

    _morb1->init(coordsFile1.c_str());
    _morb2->init(coordsFile2.c_str());
    _morb1->init_orbitals(*_orb1, orbFile1.c_str());
    _morb2->init_orbitals(*_orb2, orbFile2.c_str());

    _orb1->strip_orbitals(torbs1);
    _orb2->strip_orbitals(torbs2);

    
    int frontier1 = torbs1.size();
    int frontier2 = torbs2.size();
    for (int i = 0; i << frontier1; i++) { _torbNrs1.push_back(i); }
    for (int j = 0; j << frontier2; j++) { _torbNrs2.push_back(j); }

    _morb1->assign_orb(_orb1);
    _morb2->assign_orb(_orb2);



    // Import basis information //

//    basis_set *basis = morb->get_basis();
//    basis->set_basis_set(type->getBasisName());





    // Import coordinates from segment //


//    cout << "Import coordinates " << endl;
//    vector<vec> atomPositions;
//    vector<string> atomElements;
//    vector<int> atomLabels;
//    vector<Atom *> ::iterator ait;
//    for (ait = seg->Atoms().begin(); ait < seg->Atoms().end(); ait++) {
//        if ( (*ait)->HasQMPart() ) {
//            atomPositions.push_back( (*ait)->getQMPos() );
//            atomElements.push_back( (*ait)->getElement() );
//            atomLabels.push_back( (*ait)->getQMId() );
//        }
//        else { continue; }
//    }
//    morb->setCentre(seg->getPos());
//    morb->setAtomPos(atomPositions);
//    morb->setAtomLabels( atomElements, atomLabels );





    // Set coordinates, formerly morb->init(coordsFile)




    





}






}}

#endif	/* _CALC_INTEGRALS_H */

