/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _MD2QM_OBSERVER_H
#define	_MD2QM_OBSERVER_H

#include <votca/csg/cgobserver.h>
#include <votca/ctp/qmtopology.h>
#include <votca/tools/property.h>
#include <votca/ctp/statesaversqlite.h>
#include <votca/moo/units.h>

using namespace votca::ctp;

class MD2QMObserver
    : public CGObserver
{
public:
    MD2QMObserver();
    ~MD2QMObserver();


    void Initialize(QMTopology &qmtop, Property &opts);

    /// begin coarse graining a trajectory
    void BeginCG(Topology *top, Topology *top_atom);

    /// end coarse graining a trajectory
    void EndCG();

    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom = 0);

    void setCutoff(const double & cutoff){
        _cutoff = cutoff;
    }
    void setOut(const string &  out){
        _out=out;
    }
    //void setNNnames(string  nnnames);

    void print_nbs_to_file(QMNBList &nblist);

protected:
    QMTopology *_qmtop;
    /// nearest neighbor cut-off radius
    double _cutoff;
    ///  output streams for velocity averaging & diffusion
    ofstream _out_cont;
    ofstream _out_diff;
    string _out;
    StateSaverSQLite _save;
};

#endif	/* _MD2QM_OBSERVER_H */

