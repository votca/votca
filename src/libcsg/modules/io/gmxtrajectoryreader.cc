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

#include <cstdlib>
#include <iostream>
#include <votca/csg/topology.h>
#include "gmxtrajectoryreader.h"

namespace votca { namespace csg {

using namespace std;

bool GMXTrajectoryReader::Open(const string &file)
{
       _filename = file;
       return true;
}

void GMXTrajectoryReader::Close()
{
    close_trx(_gmx_status);
}

bool GMXTrajectoryReader::FirstFrame(Topology &conf)
{
    set_program_name("VOTCA");

#if GMX == 50

    output_env_t oenv;
    // _snew("oenv", oenv, 1);
    oenv = (output_env_t)malloc(sizeof(*oenv));
    output_env_init_default (oenv);

    if(!read_first_frame(oenv, &_gmx_status,(char*)_filename.c_str(),&_gmx_frame,TRX_READ_X | TRX_READ_V | TRX_READ_F))
        throw std::runtime_error(string("cannot open ") + _filename);
    //sfree(oenv);
    free(oenv);
#elif GMX == 45
    output_env_t oenv;
    // _snew("oenv", oenv, 1);
    oenv = (output_env_t)malloc(sizeof(*oenv));
    output_env_init_default (oenv);

    if(!read_first_frame(oenv, &_gmx_status,(char*)_filename.c_str(),&_gmx_frame,TRX_READ_X | TRX_READ_V | TRX_READ_F))
        throw std::runtime_error(string("cannot open ") + _filename);
    //sfree(oenv);
    free(oenv);
#elif GMX == 40
    if(!read_first_frame(&_gmx_status,(char*)_filename.c_str(),&_gmx_frame,TRX_READ_X  | TRX_READ_V | TRX_READ_F))
        throw std::runtime_error(string("cannot open ") + _filename);
#else
#error Unsupported GMX version
#endif

    matrix m;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            m[i][j] = _gmx_frame.box[j][i];
    conf.setBox(m);
    conf.setTime(_gmx_frame.time);
    conf.setStep(_gmx_frame.step);
    cout << endl;
    
    if(_gmx_frame.natoms != (int)conf.Beads().size())
        throw std::runtime_error("number of beads in trajectory do not match topology");

    //conf.HasPos(true);
    //conf.HasF(_gmx_frame.bF);
    
    for(int i=0; i<_gmx_frame.natoms; i++) {
        double r[3] = { _gmx_frame.x[i][XX],  _gmx_frame.x[i][YY], _gmx_frame.x[i][ZZ] };
        conf.getBead(i)->setPos(r);
        if(_gmx_frame.bF) {
            double f[3] = { _gmx_frame.f[i][XX],  _gmx_frame.f[i][YY], _gmx_frame.f[i][ZZ] };        
            conf.getBead(i)->setF(f);
        }
        if(_gmx_frame.bV) {
            double v[3] = { _gmx_frame.v[i][XX],  _gmx_frame.v[i][YY], _gmx_frame.v[i][ZZ] };
            conf.getBead(i)->setVel(v);
        }
    }
    return true;
}

bool GMXTrajectoryReader::NextFrame(Topology &conf)
{
#if GMX == 50
    output_env_t oenv;
    //_snew("oenv", oenv, 1);
    oenv = (output_env_t)malloc(sizeof(*oenv));
    output_env_init_default (oenv);
    if(!read_next_frame(oenv, _gmx_status,&_gmx_frame))
        return false;
    //sfree(oenv);
    free(oenv);
#elif GMX == 45
    output_env_t oenv;
    //_snew("oenv", oenv, 1);
    oenv = (output_env_t)malloc(sizeof(*oenv));
    output_env_init_default (oenv);
    if(!read_next_frame(oenv, _gmx_status,&_gmx_frame))
        return false;
    //sfree(oenv);
    free(oenv);
#elif GMX == 40
    if(!read_next_frame(_gmx_status,&_gmx_frame))
        return false;
#else
#error Unsupported GMX version
#endif

    matrix m;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            m[i][j] = _gmx_frame.box[j][i];
    conf.setTime(_gmx_frame.time);
    conf.setStep(_gmx_frame.step);
    conf.setBox(m);
    
    //conf.HasF(_gmx_frame.bF);
    
    for(int i=0; i<_gmx_frame.natoms; i++) {
        double r[3] = { _gmx_frame.x[i][XX],  _gmx_frame.x[i][YY], _gmx_frame.x[i][ZZ] };
        conf.getBead(i)->setPos(r);
        if(_gmx_frame.bF) {
            double f[3] = { _gmx_frame.f[i][XX],  _gmx_frame.f[i][YY], _gmx_frame.f[i][ZZ] };        
            conf.getBead(i)->setF(f);
        }
        if(_gmx_frame.bV) {
            double v[3] = { _gmx_frame.v[i][XX],  _gmx_frame.v[i][YY], _gmx_frame.v[i][ZZ] };
            conf.getBead(i)->setVel(v);
        }


    }
    return true;
}

}}
