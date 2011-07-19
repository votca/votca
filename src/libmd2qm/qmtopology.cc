#include "qmtopology.h"
#include "qmnblist.h"

QMTopology::QMTopology()
    : _db_id(0)
{}

QMTopology::~QMTopology()
{}


void QMTopology::Cleanup()
{
    vector<QMCrgUnit *>::iterator iter;
    for(iter=_crgunits.begin(); iter!=_crgunits.end(); ++iter)
        delete *iter;
    _crgunits.clear();
    _mcharges.clear();
    Topology::Cleanup();
    _db_id=0;
}

void QMTopology::Initialize(Topology& cg_top)
{
    CopyTopologyData(&cg_top);
    this->InitChargeUnits();
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

            QMCrgUnit * acrg = GetCrgUnitByName(molandtype);
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

void QMTopology::Update(Topology& cg_top)
{
    BeadContainer::iterator iter;
    BeadContainer::iterator iter_cg;

    assert(cg_top.Beads().size() == _beads.size());

    setBox(cg_top.getBox());
    _time = cg_top.getTime();
    _step = cg_top.getStep();

    iter_cg = cg_top.Beads().begin();
    for(iter=_beads.begin(); iter!=_beads.end(); ++iter, ++iter_cg) {
        (*iter)->setPos((*iter_cg)->getPos());
        if((*iter_cg)->HasU())
            (*iter)->setU((*iter_cg)->getU());
        if((*iter_cg)->HasV())
            (*iter)->setV((*iter_cg)->getV());
        if((*iter_cg)->HasW())
            (*iter)->setW((*iter_cg)->getW());
        QMBead * b = (QMBead*)(*iter);

        if((*iter_cg)->getSymmetry() == 1) {
            (*iter)->setU(vec(1,0,0));
            (*iter)->setV(vec(0,1,0));
            (*iter)->setW(vec(0,0,1));
        }

        b->UpdateCrg();
    }
}

void QMTopology::LoadListCharges(const string &file)
{
    _jcalc.Initialize(file);
}


Molecule *QMTopology::AddAtomisticBeads(CrgUnit *crg, Topology * totop){
     
    mol_and_orb * atoms = crg->rotate_translate_beads();
    totop->CreateResidue("DUM");
    Molecule *mi = totop->CreateMolecule((crg)->getName());
    mi->setUserData(crg);   //mi->getUserData<CrgUnit>();
    for (int i=0;i<atoms->getN();i++){
        vec pos = unit<bohr,nm>::to(atoms->GetPos(i));
        string atomtype = string( atoms->gettype(i) ); //+ string("-") + lexical_cast<string>(crg->getId());
        BeadType * bt= totop->GetOrCreateBeadType(atomtype);
        Bead * bead = totop ->CreateBead(1, atomtype,bt,0, 0, 0.);
        bead->setPos(pos);
        mi->AddBead(bead,"???");
        if(crg->getType()->GetCrgUnit().getChargesNeutr()) {
            double charge_of_bead_neutral=crg->getType()->GetCrgUnit().getChargesNeutr()->mpls[i];
            bead->setQ(charge_of_bead_neutral);
        }
      
    }
    delete atoms->getorbs();
    delete atoms;
    return mi;
}

Bead *QMTopology::CreateBead(byte_t symmetry, string name, BeadType *type, int resnr, double m, double q)
{
    QMBead *bead = new QMBead(this, _beads.size(), type, symmetry, name, resnr, m, q);
    _beads.push_back(bead);
    return bead;
}

// TODO: this function should not be in qmtopology!
void QMTopology::ComputeAllTransferIntegrals(){
    for(QMNBList::iterator iter = _nblist.begin();
        iter!=_nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->Crg1();
        CrgUnit *crg2 = (*iter)->Crg2();
        vector <double> Js = GetJCalc().CalcJ(*crg1, *crg2);
        (*iter)->setJs(Js);
    }
}

   // In topology.cc???
//Copy charges to either charged or neutral case

void QMTopology::CopyCharges(CrgUnit *crg, Molecule *mol)
{
    if(mol->BeadCount() != crg->getType()->GetCrgUnit().getN())
        throw std::runtime_error("QMTopology::CopyCharges: number of atoms in crgunit does not match number of beads in molecule");

    if(!crg->getType()->GetCrgUnit().getChargesNeutr())
        throw std::runtime_error("QMTopology::CopyCharges: no charges defined");

    //loop over all beads in that molecule
    for (int i = 0; i < mol->BeadCount(); i++) {
        //get charge
        double charge_of_bead_i_neutral = crg->getType()->GetCrgUnit().getChargesNeutr()->mpls[i];
        //set charge
        mol->getBead(i)->setQ(charge_of_bead_i_neutral);
    }
}

//Copy charges to either charged or neutral case
void QMTopology::CopyChargesOccupied(CrgUnit *crg, Molecule *mol)
{
    if(mol->BeadCount() != crg->getType()->GetCrgUnit().getN())
        throw std::runtime_error("QMTopology::CopyCharges: number of atoms in crgunit does not match number of beads in molecule");

    if(!crg->getType()->GetCrgUnit().getChargesCrged())
        throw std::runtime_error("QMTopology::CopyCharges: no charges defined");

    
    if(crg->getType()->GetCrgUnit().getChargesCrged()->mpls.size() == 0)
        throw std::runtime_error("QMTopology::CopyCharges: no charges defined (mpls.size()==0, that is strange!)");

    //loop over all beads in that molecule
    for (int i = 0; i < mol->BeadCount(); i++) {
        //get charge
        double charge_of_bead_i_charged = crg->getType()->GetCrgUnit().getChargesCrged()->mpls[i];
        //set charge
        mol->getBead(i)->setQ(charge_of_bead_i_charged);
    }
}

QMCrgUnit *QMTopology::CreateCrgUnit(const string &name, const string &type_name, int molid)
{
    CreateCrgUnit(_crgunits_by_id.size()+1, name, type_name, molid);
}

QMCrgUnit *QMTopology::CreateCrgUnit(int id, const string &name, const string &type_name, int molid)
{
    map<int, QMCrgUnit*>::iterator iter;
    iter = _crgunits_by_id.find(id);
    if(iter != _crgunits_by_id.end())
        throw std::runtime_error("charge unit with id " + lexical_cast<string>(id) + " already exists");

    if(GetCrgUnitByName(name))
        throw std::runtime_error("charge unit with name " + name + " already exists");

    QMCrgUnit *crg;

    CrgUnitType *type = _jcalc.GetCrgUnitTypeByName(type_name);
    if(!type)
        throw runtime_error("Charge unit type not found: " + type_name);

    crg = new QMCrgUnit(id, type, molid);

    _mcharges.insert(make_pair(name, crg));
    _crgunits.push_back(crg);
    _crgunits_by_id[id] = crg;
    crg->setName(name);
    return crg;
}
