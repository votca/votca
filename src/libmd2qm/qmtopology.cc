#include "qmtopology.h"
#include "qmnblist.h"

QMTopology::QMTopology()
{}

QMTopology::~QMTopology()
{}


void QMTopology::Initialize(Topology& cg_top)
{
    CopyTopologyData(&cg_top);
    this->InitChargeUnits();
}

void QMTopology::Update(Topology& cg_top)
{
    BeadContainer::iterator iter;
    BeadContainer::iterator iter_cg;

    assert(cg_top.Beads().size() == _beads.size());

    _box = cg_top.getBox();
    _time = cg_top.getTime();
    _step = cg_top.getStep();

    iter_cg = cg_top.Beads().begin();
    for(iter=_beads.begin(); iter!=_beads.end(); ++iter, ++iter_cg) {
        (*iter)->setPos((*iter_cg)->getPos());
        (*iter)->setU((*iter_cg)->getU());
        (*iter)->setV((*iter_cg)->getV());
        (*iter)->setW((*iter_cg)->getW());
        QMBead * b = (QMBead*)(*iter);
        b->UpdateCrg();
    }
}

void QMTopology::LoadListCharges(const string &file)
{
    _jcalc.Init(file);
}


void QMTopology::AddAtomisticBeads(CrgUnit *crg, Topology * totop){
    
    mol_and_orb * atoms = crg->rotate_translate_beads();
    totop->CreateResidue("DUM");
    for (int i=0;i<atoms->getN();i++){
        vec pos = atoms->GetPos(i) *(RA/10.) ;//convert to nm!
        string atomtype = string( atoms->gettype(i) )+ string("-") + lexical_cast<string>(crg->GetId());
        BeadType * bt= totop->GetOrCreateBeadType(atomtype);
        Bead * bead = totop ->CreateBead(1, atomtype,bt,0, 0, 0.);
        bead->setPos(pos); 
    }
    delete atoms;
}

Bead *QMTopology::CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q)
{
    QMBead *bead = new QMBead(this, _beads.size(), type, symmetry, name, resnr, m, q);
    _beads.push_back(bead);
    return bead;
}

void QMTopology::InitChargeUnits(){
    BeadContainer::iterator itb;
    for (itb = _beads.begin() ; itb< _beads.end(); ++itb){
        QMBead * bead = dynamic_cast<QMBead *>(*itb);
        //initialise the crgunit * only if appropriate extra info is in the cg.xml file
        if ( (bead->Options()).exists("qm.crgunitname")){
            string namecrgunittype = bead->getType()->getName();
            int intpos = (bead->Options()).get("qm.position").as<int>();
            string namecrgunit = (bead->Options()).get("qm.crgunitname").as<string>();

            CrgUnitType *crgtype  = _jcalc.GetCrgUnitTypeByName(namecrgunittype);

            //determine whether it  has been created already
            int molid= bead->getMolecule()->getId();
            string molandtype = lexical_cast<string>(molid)+":"+namecrgunit;
            map <string, CrgUnit*>::iterator  itm= _mcharges.find(molandtype);
            if (itm == _mcharges.end()){

                CrgUnit * acrg = new CrgUnit(_lcharges.size(), crgtype, molid);
                _mcharges.insert(make_pair(molandtype, acrg));
                _lcharges.push_back(acrg);
                bead->setCrg(acrg);
                bead->setiPos(intpos);
            }
        }
        else{
            bead->setCrg(NULL);
        }
    }
}

void QMTopology::ComputeAllTransferIntegrals(){
    for(QMNBList::iterator iter = _nblist.begin();
        iter!=_nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
        vector <double> Js = GetJCalc().GetJ(*crg1, *crg2);
        (*iter)->setJ(Js);
    }
}