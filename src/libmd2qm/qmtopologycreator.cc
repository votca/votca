#include "qmtopologycreator.h"


void QMTopologyCreator::Initialize(Topology& cg_top)
{
    CopyTopologyData(&cg_top);
    this->InitChargeUnits();
}

void QMTopologyCreator::InitChargeUnits(){
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
                acrg = CreateCrgUnit(0, molandtype, namecrgunittype, molid);

            bead->setCrg(acrg);
            bead->setiPos(intpos);
        }
        else{
            bead->setCrg(NULL);
        }
    }
}