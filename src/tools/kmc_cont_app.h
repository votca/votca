/* 
 * File:   kmc_cont_app.h
 * Author: vehoff
 *
 * Created on March 4, 2010, 4:50 PM
 */

#ifndef _KMC_CONT_APP_H
#define	_KMC_CONT_APP_H

#include "qmapplication.h"
#include <kmc/vertex.h>
#include <kmc/hoppers.h>
#include <kmc/kmc.h>
#include <kmc/graph.h>

class KmcCont : public QMApplication
{
public:
    KmcCont();
    ~KmcCont();

    void HelpText();
    void Initialize();
    bool EvaluateFrame();

private:
    /// electric field
    vec _E;
    /// total simulation time
    double _total_time;
    /// time step
    double _dt;
    /// number of charges simultaneously present in a simultion
    int _ncrg;
    /// number of KMC runs to be performed
    int _nruns;
    ///  output streams for velocity averaging & diffusion
    ofstream _out_cont;
    ofstream _out_diff;

    /// creation of KMC graph

    void make_kmc_graph(graph *a, QMNBList &nblist);
};

#endif	/* _KMC_CONT_APP_H */

