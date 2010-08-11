/* 
 * File:   ctp_connectivity.cc
 * Author: schrader
 *
 * Created on August 10, 2010, 12:24 PM
 */

#include <stdlib.h>
#include "qmapplication.h"
#include <calculatorfactory.h>

using namespace std;

/*
 * Connectivity Graph
 */
class QMAppConnectivity : public QMApplication
{
public:
    void HelpText() {}

    void AddSpecificOptions(){
        _op_desc_specific.add_options()
        ("min", boost::program_options::value<double>(), "Minimum value")
        ("max", boost::program_options::value<double>(), "Maximum value")
        ("scale1", boost::program_options::value<double>(), "Drawing scale factor for bond thickness")
        ("scale2", boost::program_options::value<double>(), "Drawing scale factor for spheres")
        ;
    }

    void BeginEvaluate(){
        //TODO: check if options are actually set, read in from xml file
        label=0;
        min=_op_vm["min"].as<double>();
        max=_op_vm["max"].as<double>();
        scale1=_op_vm["scale1"].as<double>();
        scale2=_op_vm["scale2"].as<double>();
    }

    bool EvaluateFrame(){

        //TODO: make numbering according to frame number
        label++;

        //TODO: cleanup following code
        FILE* outfile;
        FILE* outvmd;
        char namepdb[50];
        char namevmd[50];
        sprintf(namepdb, "connectivity_%04d.pdb", label);
        sprintf(namevmd, "connectivity_%04d.vmd", label);
        outfile = fopen(namepdb,"w");
        outvmd = fopen(namevmd,"w");

        string type;
        map < pair <int,int>, double >::iterator it_map; // map of interesting nbs
        pair <int,int> poi; // pair of interest
        double Joi; // transfer integral of interest
        map <pair <int,int>, double > nbmap; // neighbours of interest
        CrgUnit* nb1; // 1st charge unit
        CrgUnit* nb2; // 2nd charge unit
        vec com; // center of mass
        double dist; // distance between molecules ignoring PBCs

        fprintf(outfile, "CONNECTIVITY GRAPH \n");
        type = "C"; // using Carbon atoms for the centers of mass

        // CrgUnits to pdb file
        vector<CrgUnit *> lcharges = _qmtop.CrgUnits();
        for (vector<CrgUnit *>::iterator it_crgs = lcharges.begin(); it_crgs != lcharges.end(); it_crgs++){
                com = (*it_crgs)->GetCom();
                //if(zllim<com.getZ() && com.getZ()<zulim && xllim<com.getX() && com.getX()<xulim && yllim<com.getY() && com.getY()<yulim)
                    fprintf(outfile, "ATOM  %5d  %3s %3s %5d     %7.3f %7.3f %7.3f\n", (*it_crgs)->getId(), type.c_str(), "CG", (*it_crgs)->getId(), com.getX(), com.getY(), com.getZ());
        }

        // Transfer integrals in given range
        QMNBList &nblist = _qmtop.nblist();
        for(QMNBList::iterator it_neigh = nblist.begin(); it_neigh!=nblist.end(); it_neigh++){
            if( log10((*it_neigh)->calcJeff2()) > min && log10((*it_neigh)->calcJeff2()) < max ){
                nb1 = _qmtop.GetCrgUnit((*it_neigh)->Crg1()->getId());
                nb2 = _qmtop.GetCrgUnit((*it_neigh)->Crg2()->getId());
                dist = abs(nb1->GetCom()-nb2->GetCom());
                if (dist <= (*it_neigh)->dist()) //compare naive distance with distance regarding pbc to ignore bonds over pbc
                {
                //if (dist < _d_nn && zllim<nb1->GetCom().getZ() && nb1->GetCom().getZ()<zulim && xllim<nb1->GetCom().getX()  && nb1->GetCom().getX()<xulim  && yllim<nb1->GetCom().getY()  && nb1->GetCom().getY()<yulim){
                    poi.first = nb1->getId();
                    poi.second = nb2->getId();
                    //use logarithmic thickness scaling
                    Joi = -1/log10((*it_neigh)->calcJeff2());// sqrt( (*it_neigh)->calcJeff2() );
                    nbmap.insert( pair < pair<int,int>, double > (poi, Joi));
                    fprintf(outfile, "CONECT %4d %4d \n", nb1->getId(), nb2->getId());
                //}
                }
            }
        }
        fclose(outfile);

        // vmd file
        fprintf(outvmd, "display projection Orthographic\n\n");
        fprintf(outvmd, "mol new %s type pdb waitfor all\n", namepdb);
        fprintf(outvmd, "color Display {Background} white\n");
        fprintf(outvmd, "trans 100 00 0\n");
        fprintf(outvmd, "axes location off\n\n");
        fprintf(outvmd, "# first define all useful selections\n\n");

        fprintf(outvmd, "set all \"{all}\"\n");

        int i=0;
        for ( it_map = nbmap.begin(); it_map != nbmap.end(); ++it_map ){
            fprintf(outvmd, "set nbs_%d \"{resid %d or resid %d}\"\n", i, (*it_map).first.first, (*it_map).first.second);
            i++;
        }

        fprintf(outvmd, "\n\n");
        fprintf(outvmd, "set all_rep \"VDW %3.2f 10\"\n", scale2);

        i=0;
        for ( it_map = nbmap.begin(); it_map != nbmap.end(); ++it_map ){
            fprintf(outvmd, "set nbs_%d_rep \"bonds %3.2f 10\"\n", i, fabs((*it_map).second)*scale1);
            i++;
        }

        fprintf(outvmd, "\n");
        fprintf(outvmd, "# delete old representation\n\n");
        fprintf(outvmd, "mol delrep 0 top\n\n");

        fprintf(outvmd, "#create new representation\n\n");
        fprintf(outvmd, "mol representation $all_rep\n");
        fprintf(outvmd, "mol color ColorID 2\n");
        fprintf(outvmd, "mol selection $all\n");
        fprintf(outvmd, "mol material Opaque\n");
        fprintf(outvmd, "mol addrep top\n\n");

        i=0;
        for ( it_map = nbmap.begin(); it_map != nbmap.end(); ++it_map ){
            fprintf(outvmd, "mol representation $nbs_%d_rep \n", i);
            if((*it_map).second<0.0){fprintf(outvmd, "mol color ColorID 1\n");}
            else{fprintf(outvmd, "mol color ColorID 0\n");}
            fprintf(outvmd, "mol selection $nbs_%d\n", i);
            fprintf(outvmd, "mol material Opaque\n");
            fprintf(outvmd, "mol addrep top\n\n");
            i++;
        }

        fclose(outvmd);
    }

    void EndEvaluate(){
    }

    private:
    int label;
    double min, max;
    double scale1,scale2;
};

int main(int argc, char** argv) {
    QMAppConnectivity qmappcon;
    qmappcon.Run(argc, argv);
    return (EXIT_SUCCESS);
}

