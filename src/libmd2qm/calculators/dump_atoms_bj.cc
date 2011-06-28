#include <stdlib.h>
#include "dump_atoms_bj.h"
#include <math.h>
#include <list>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

void DumpAtomsBJ::Initialize(QMTopology *top, Property *options) {
    _options = options;
    if (options->exists("options.dump_atoms_params.cutoff")) {
        _dump_cutoff = options->get("options.dump_atoms_params.cutoff").as<double>();
        cout << "Writing out atomic XYZ coordinates for molecular cutoff: " << _dump_cutoff << " nm" << endl;
    } else {
        _dump_cutoff = 50.0;
        cout << "Warning: No cutoff for molecules has been provided, using default cutoff of 50nm." << endl;
    }

}

bool DumpAtomsBJ::EvaluateFrame(QMTopology *top) {
    QMNBList& nblist = top->nblist();
    for (QMNBList::iterator ipair = nblist.begin(); ipair != nblist.end(); ++ipair) {

        QMPair *pair = *ipair;
        QMCrgUnit *crg1 = pair->Crg1();
        QMCrgUnit *crg2 = pair->Crg2();

        vector<QMCrgUnit *> lcharges = top->CrgUnits();
        Topology atop;
        atop.setBox(top->getBox());
        //Add the beads of charge units i and j to the atop
        Molecule *molcrg1 = top->AddAtomisticBeads(crg1, &atop);
        Molecule *molcrg2 = top->AddAtomisticBeads(crg2, &atop);
        //Add the rest of the beads to atop
        for (vector<QMCrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++) {
            if (crg1->getId() == (*itl)->getId()) continue;
            if (crg2->getId() == (*itl)->getId()) continue;
            top->AddAtomisticBeads(*itl, &atop);
        }

        //vec diff = 0.5 * (crg1->GetCom() + crg2->GetCom());
        //vec center = crg2->GetCom() + diff;
        vec center = 0.5 * (crg1->GetCom() + crg2->GetCom());
        //write out atoms of molecule A=i
        string filename = "pair_" + boost::lexical_cast<string > (crg1->getId() + 1) + "_" + boost::lexical_cast<string > (crg2->getId() + 1) + "_A";
        FILE * data;
        data = fopen(filename.c_str(), "w");
        //Loop over all beads in crgunit A=i

         MoleculeContainer::iterator imol;
         for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
          if (crg->getId() != crg1->getId()) continue;
          vec bcs = atop.BCShortestConnection(center,crg->GetCom());
          vec dist = crg->GetCom() - center ;
          vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }
        //loop over all others except B=j and A=i
         for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
            vec bcs = atop.BCShortestConnection(center,crg->GetCom());
            if (abs(bcs) > _dump_cutoff) continue;
            if (crg->getId() == crg1->getId()) continue;
            if (crg->getId() == crg2->getId()) continue;
                 vec dist = crg->GetCom() - center ;
                 vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }
        fclose(data);


        //write out all atoms of molecule B=j
        string filename2 = "pair_" + boost::lexical_cast<string > (crg1->getId() + 1) + "_" + boost::lexical_cast<string > (crg2->getId() + 1) + "_B";
                FILE * data2;
        data2 = fopen(filename2.c_str(), "w");
        //Loop over all beads in crgunit B=j
         for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
          if (crg->getId() != crg2->getId()) continue;
          vec bcs = atop.BCShortestConnection(center,crg->GetCom());
          vec dist = crg->GetCom() - center ;
          vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data2, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }
        //loop over all others except B=j and A=i
         for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
            vec bcs = atop.BCShortestConnection(center,crg->GetCom());
            if (abs(bcs) > _dump_cutoff) continue;
            if (crg->getId() == crg1->getId()) continue;
            if (crg->getId() == crg2->getId()) continue;
                 vec dist = crg->GetCom() - center ;
                 vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data2, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }
        fclose(data2);

            //write out all atoms of molecule A=i and molecule B=j
        string filename3 = "pair_" + boost::lexical_cast<string > (crg1->getId() + 1) + "_" + boost::lexical_cast<string > (crg2->getId() + 1) + "_AB";
                FILE * data3;
        data3 = fopen(filename3.c_str(), "w");
        //Loop over all beads in crgunit A=i
        for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
          if (crg->getId() != crg1->getId()) continue;
          vec bcs = atop.BCShortestConnection(center,crg->GetCom());
          vec dist = crg->GetCom() - center ;
          vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data3, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }
         //Loop over all beads in crgunit B=j
         for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
          if (crg->getId() != crg2->getId()) continue;
          vec bcs = atop.BCShortestConnection(center,crg->GetCom());
          vec dist = crg->GetCom() - center ;
          vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data3, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }

        //loop over all others except A=i and B=j
         for (imol = atop.Molecules().begin(); imol != atop.Molecules().end(); imol++) {
          Molecule *mol = *imol;
          CrgUnit *crg = mol->getUserData<CrgUnit > ();
          vec bcs = atop.BCShortestConnection(center,crg->GetCom());
            if (abs(bcs) > _dump_cutoff) continue;
            if (crg->getId() == crg1->getId()) continue;
            if (crg->getId() == crg2->getId()) continue;
                
          vec dist = crg->GetCom() - center ;
          vec diffs = bcs - dist;
               for (int j = 0; j != mol->BeadCount(); j++) {
                Bead *bj = mol->getBead(j);
                vec r_v = (bj->getPos())+diffs-center;
                fprintf(data3, "%.6f %.6f %.6f\n", r_v.getX()*10.0, r_v.getY()*10.0, r_v.getZ()*10.0);
            }
        }
        fclose(data3);
    }
    
    return true;
}

   


