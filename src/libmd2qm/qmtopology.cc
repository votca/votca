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

void QMTopology::Cleanup()
{
    list<CrgUnit *>::iterator iter;
    for(iter=_lcharges.begin(); iter!=_lcharges.end(); ++iter)
        delete *iter;
    _lcharges.clear();
    _mcharges.clear();
    Topology::Cleanup();
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
    _jcalc.Initialize(file);
}


void QMTopology::AddAtomisticBeads(CrgUnit *crg, Topology * totop){
    
    mol_and_orb * atoms = crg->rotate_translate_beads();
    totop->CreateResidue("DUM");
    for (int i=0;i<atoms->getN();i++){
        vec pos = unit<bohr,nm>::to(atoms->GetPos(i));
        string atomtype = string( atoms->gettype(i) )+ string("-") + lexical_cast<string>(crg->getId());
        BeadType * bt= totop->GetOrCreateBeadType(atomtype);
        Bead * bead = totop ->CreateBead(1, atomtype,bt,0, 0, 0.);
        bead->setPos(pos); 
    }
    delete atoms->getorbs();
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

            //determine whether it  has been created already
            int molid= bead->getMolecule()->getId();
            string molandtype = lexical_cast<string>(molid)+":"+namecrgunit;
            
            CrgUnit * acrg = GetCrgUnitByName(molandtype);
            if(acrg == NULL)
                acrg = CreateCrgUnit(molandtype, namecrgunittype, molid);
            
            bead->setCrg(acrg);
            bead->setiPos(intpos);
        }
        else{
            bead->setCrg(NULL);
        }
    }
}

void QMTopology::ComputeAllTransferIntegrals(){
    for(QMNBList::iterator iter = _nblist.begin();
        iter!=_nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->Crg1();
        CrgUnit *crg2 = (*iter)->Crg2();
        vector <double> Js = GetJCalc().CalcJ(*crg1, *crg2);
        (*iter)->setJs(Js);
    }
}

void QMTopology::ComputeAllElectrostaticEnergies(const double& epsilon){
    list <CrgUnit * >::iterator itcrg = _lcharges.begin();
    list <CrgUnit * >::iterator itneutr= _lcharges.begin();
    for ( itcrg = _lcharges.begin() ; itcrg != _lcharges.end() ; ++itcrg){
        double nrg=0.;
        for (itneutr = _lcharges.begin() ; itneutr != _lcharges.end() ; ++itneutr){
            if ( itcrg != itneutr){
//              nrg += ENERGY FOR ITCRG WITH ITNEUTR WITH ITCRG CHARGED;
//              NRG -= ENERGY FOR ITCRG WITH ITNEUTR WITH ITCRG NEUTRAL;

                QMPair twocharges(*itcrg, *itneutr, this);
                CrgUnit * crg1= twocharges.first;
                CrgUnit * crg2= twocharges.second;

                double contr = _jcalc.EstaticDifference(*crg1, *crg2);
                nrg += contr;

      //          clog << "[kmcdata.cc] contribution from" << (*itneutr)->getType()->GetName() << " " << (*itneutr)->getId() <<
      //                 " to " << (*itcrg)->getType()->GetName() << " " << (*itcrg)->getId() <<
      //                 " distance " << abs(r_ij) << " nrg " << contr << endl;


            }
        }
        nrg /= epsilon;
        clog << "[kmc_data.cc] FINALE ENERGY for "
                << (*itcrg)->GetCom() << " " << (*itcrg)->getId() << " " << (*itcrg)->getType()->GetName()
                << " " << nrg << endl;
        (*itcrg)->setEnergy((*itcrg)->getEnergy()+ nrg);
    }
}
