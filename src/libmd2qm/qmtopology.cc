#include "qmtopology.h"
#include "qmnblist.h"

QMTopology::QMTopology()
{
    _nblist = NULL;
}

QMTopology::~QMTopology()
{
    if(_nblist)
        delete _nblist;
    _nblist = NULL;
}


void QMTopology::Initialize(Topology& cg_top)
{
    CopyTopologyData(&cg_top);
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
    for(iter=_beads.begin(); iter!=_beads.end(); ++iter) {
        (*iter)->setPos((*iter_cg)->getPos());
        (*iter)->setU((*iter_cg)->getU());
        (*iter)->setV((*iter_cg)->getV());
        (*iter)->setW((*iter_cg)->getW());
        QMBead * b = dynamic_cast<QMBead*>(*iter);
        b->QMBead::UpdateCrg();
    }
}

void QMTopology::LoadListCharges(const string &file)
{
    _jcalc.Init(file);
}


void QMTopology::AddAtomisticBeads(CrgUnit *crg, Topology * totop){
    
    mol_and_orb * atoms = crg->rotate_translate_beads();

    for (int i=0;i<atoms->getN();i++){
        vec pos = atoms->GetPos(i);
        string atomtype = "QMAT-"+string( atoms->gettype(i) );
        BeadType * bt= totop->GetOrCreateBeadType(atomtype);
        int nbead = totop-> BeadCount();
        Bead * bead = totop ->CreateBead(1, atomtype,bt,nbead, 0, 0.);
        bead->setPos(pos);
    }
    delete atoms;
}

Bead *QMTopology::CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q)
{
    QMBead *bead = new QMBead(this, _beads.size(), type, symmetry, name, resnr, m, q);
    _beads.push_back(bead);

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
        if (itm != _mcharges.end()){
            
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
    return bead;
}
