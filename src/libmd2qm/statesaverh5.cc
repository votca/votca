#include <statesaverh5.h>
#include "hdf5_hl.h"

void StateSaverH5::Open(QMTopology& qmtop, const string& file)
{
    _file = H5Fcreate(file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    _frame = 0;
    _qmtop=&qmtop;
}

void StateSaverH5::Close()
{
    herr_t status;
    status = H5Fclose(_file);
}

void StateSaverH5::WriteFrame()
{
    string str="frame"+lexical_cast<string>(_frame);
    hid_t frame = H5Gcreate(_file, str.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    write(frame, "time", _qmtop->getStep());
    write(frame, "step", _qmtop->getTime());
    write(frame, "box", _qmtop->getBox());

    WriteMolecules(frame);
    WriteBeads(frame);
    WritePairs(frame);

    H5Gclose(frame);
    _frame++;
}

void StateSaverH5::WriteMolecules(hid_t context)
{
    hid_t mols = H5Gcreate(context, "molecules", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write(mols, "count", _qmtop->MoleculeCount());

    int imol=0;
    for (MoleculeContainer::iterator iter = _qmtop->Molecules().begin();
            iter != _qmtop->Molecules().end() ; ++iter) {
        Molecule *mol=*iter;
        string str=lexical_cast<string>(imol++);
        hid_t mol_grp = H5Gcreate(mols, str.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        write(mol_grp, "id", mol->getId());
        write(mol_grp, "name", mol->getName());
        H5Gclose(mol_grp);
    }
    H5Gclose(mols);
}

void StateSaverH5::WriteBeads(hid_t context) {
    hid_t beads = H5Gcreate(context, "beads", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write(beads, "count", _qmtop->BeadCount());
    int ibead=0;
    for (BeadContainer::iterator iter = _qmtop->Beads().begin();
            iter != _qmtop->Beads().end(); ++iter) {
        QMBead *bi = dynamic_cast<QMBead*> (*iter);
        string str=lexical_cast<string>(ibead++);
        hid_t bead_grp = H5Gcreate(beads, str.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        write(bead_grp, "name", bi->getName());
        write(bead_grp, "symmetry", (int)bi->getSymmetry());
        write(bead_grp, "type", bi->getType()->getName());
        write(bead_grp, "resnr", bi->getResnr());
        write(bead_grp, "mass", bi->getM());
        write(bead_grp, "charge", bi->getQ());

        write(bead_grp, "crgunit", bi->GetCrgUnit()->getName());
        write(bead_grp, "energy", bi->GetCrgUnit()->getEnergy());
        write(bead_grp, "crgunit_index", bi->getiPos());
        write(bead_grp, "pos", bi->getPos());
        write(bead_grp, "u", bi->getU());
        write(bead_grp, "v", bi->getV());
        write(bead_grp, "molid", bi->getMolecule()->getId());

        H5Gclose(bead_grp);
    }
    H5Gclose(beads);
}

void StateSaverH5::WritePairs(hid_t context) {
    hid_t pairs = H5Gcreate(context, "pairs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    write(pairs, "count", (int)_qmtop->nblist().size());
    int ipair=0;

    QMNBList &nblist = _qmtop->nblist();

    //cout <<"There are so many pairs in nblist: " <<(int)nblist.size()<<"\n";
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        string str=lexical_cast<string>(ipair++);
        hid_t pair_grp = H5Gcreate(pairs, str.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        QMPair *pair = *iter;
        QMCrgUnit *crg1 = (*iter)->first;
        QMCrgUnit *crg2 = (*iter)->second;

        write(pair_grp, "crgunitid1", (int)crg1->getId());
        write(pair_grp, "crgunitid2", (int)crg2->getId());
        write(pair_grp, "crgunit1", crg1->getName());
        write(pair_grp, "crgunit2", crg2->getName());
        //write(pair_grp, "integrals", pair->Js());
        write(pair_grp, "rate12", pair->rate12());
        write(pair_grp, "rate21", pair->rate21());
        write(pair_grp, "r", pair->r());
        H5Gclose(pair_grp);
    }
    H5Gclose(pairs);
}
