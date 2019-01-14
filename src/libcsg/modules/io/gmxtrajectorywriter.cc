/* 
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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


#include <string>
#include "gmxtrajectorywriter.h"

namespace votca { namespace csg {

using namespace std;

void GMXTrajectoryWriter::Open(string file, bool bAppend)
{
    //char c[1] = bAppend ? "a" : "w";
    _file = open_trx((char *)file.c_str(), "w");
}

void GMXTrajectoryWriter::Close()
{
    close_trx(_file);
}

void GMXTrajectoryWriter::Write(Topology *conf)
{
    static int step=0;   
    int N = conf->BeadCount();
    t_trxframe frame;
    rvec *x = new rvec[N];
    rvec *v = NULL;
    rvec *f = NULL;
    matrix box = conf->getBox();
    
    frame.natoms = N;
    frame.bTime = true;
    frame.time = conf->getTime();
    frame.bStep = true;
    frame.step = conf->getStep();;
    frame.x = x;
    frame.bLambda=false;
    frame.bAtoms=false;
    frame.bPrec=false;
    frame.bX = true;
    frame.bF=conf->HasForce();
    frame.bBox=true;
    frame.bV=conf->HasVel();

    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            frame.box[j][i] = box[i][j];
    
for(int i=0; i<N; ++i) {
        vec pos = conf->getBead(i)->getPos();
        x[i][0] = pos.getX();
        x[i][1] = pos.getY();
        x[i][2] = pos.getZ();
    }

if (frame.bV){
    v = new rvec[N];
    for(int i=0; i<N; ++i) {
        frame.v = v;
        vec vel = conf->getBead(i)->getVel();
        v[i][0] = vel.getX();
        v[i][1] = vel.getY();
        v[i][2] = vel.getZ();
    }
}
 if (frame.bF){
     f = new rvec[N];
    for(int i=0; i<N; ++i) {
        frame.f = f;
        vec force = conf->getBead(i)->getF();
        f[i][0] = force.getX();
        f[i][1] = force.getY();
        f[i][2] = force.getZ();
    }
}
     
    write_trxframe(_file, &frame, NULL);

    step++;
    delete[] x;
    if (frame.bV) delete[] v;
    if (frame.bF) delete[] f;
}

}}
