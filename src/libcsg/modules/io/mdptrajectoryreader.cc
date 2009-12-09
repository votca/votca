/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
#include "mdptrajectoryreader.h"
#include "configuration.h"

bool MDPTrajectoryReader::Open(const string &file)
{
    
    _fl = fopen(file.c_str(), "rb");
    if(!_fl)
        throw string("cannot open " + file);
}

void MDPTrajectoryReader::Close()
{
    fclose(_fl);
    _nmols.clear();
    _natoms.clear();
}


bool MDPTrajectoryReader::FirstFrame(Configuration &conf)
{
    FILE *fl = _fl;
    long I4;
    float R4;
    double R8;
        
    conf.HasPos(true);
/* 
 *   The header line looks as follows (in fortran binary out, so with record blocks):
 *   
 *   I4    0/1 has velocities?
 *   R8    dt
 *   3*R8  box size (Angstr)
 *   R8    unit in m (10^-10 for Angstr)
 *   I4    number if mol types
 *   for each molecule:
 *   I4    nspec (number of molecules of this type)
 *   I4    nsites (number of atoms in molecule)
 *   L4    true ??
 *
 * */

// read record begin
    fread(&I4, sizeof(I4), 1, fl);
    
    // read ivel
    fread(&I4, sizeof(I4), 1, fl);
    
    double dt, unit;
    vec box;
    
    fread(&dt, sizeof(dt), 1, fl);
    fread(&box, sizeof(box), 1, fl);
    fread(&unit, sizeof(unit), 1, fl);
    fread(&_moltypes, sizeof(_moltypes), 1, fl);
    
            
    for(int i=0; i<_moltypes; ++i) {
        fread(&I4, sizeof(I4), 1, fl);
        _nmols.push_back(I4);
        fread(&I4, 1, sizeof(I4), fl);
        _natoms.push_back(I4);
        fread(&I4, sizeof(I4), 1, fl);
    }

    // read record end
    fread(&I4, sizeof(I4), 1, fl);
    // now read in the molecule types        
    /*
     *  for each type 
     *    for each atom
     *      R8       mass
     *  I4 ?         List (0: has coordinates, 1: no coordinates)
     */
      
    fread(&I4, sizeof(I4), 1, fl);
     
    for(int type=0; type<_moltypes; ++type) {
        for(int atom=0; atom<_natoms[type]; ++atom) {
            fread(&R8, sizeof(R8), 1, fl);
            cout << R8 << endl;
        }

        // is molecule read in trajectory?
        fread(&I4, sizeof(I4), 1, fl);
        cout << I4 << endl;
    }
    
    // read record end
    fread(&I4, sizeof(I4), 1, fl);
    
    return NextFrame(conf);
}

bool MDPTrajectoryReader::NextFrame(Configuration &conf)
{
    long I4;
    float R4;
    double R8;
    FILE *fl = _fl;
    
    if(feof(fl)) return false;
    
    /* record
     *
     *  I4        ivel, 0/1, has velocity?
     *  R8        time
     *  R8        comtime (averages, ..)?
     *  R8        temperature
     *  R8        pressure
     *  R8        Epot
     *  3*R8      box size
     *  L4        list: included in file?
     *
     */

    long iVel;
    
    long rec;
    // read in record marker, to skip list, was buggy in fotran code...
    
    fread(&rec, sizeof(rec), 1, fl);  // begin record
    fread(&iVel, sizeof(iVel), 1, fl);
    conf.HasVel(iVel == 1 ? true : false);
    
    double t;
    fread(&t, sizeof(t), 1, fl);
    conf.setTime(t);
    
    fread(&R8, sizeof(R8), 1, fl);  // comtime
    fread(&R8, sizeof(R8), 1, fl);  // temperature
    fread(&R8, sizeof(R8), 1, fl);  // pressure
    fread(&R8, sizeof(R8), 1, fl);  // Epot
    
    vec box;
    fread(&box, sizeof(box), 1, fl);
    cout << box << endl;
    
    conf.setBox(matrix(vec(box.x(), 0, 0), vec(0, box.y(), 0), vec(0, 0, box.y())));
    I4 = rec + 1;
    while(I4!=rec)  // skip everything to end of record
        fread(&I4, sizeof(I4), 1, fl);  // end record

    
    /* for each moleculetype record (if list)
     *  record
     *  3*R4*natoms*ntypes      pos
     *
     *  record if ivel
     *  3*R4*natoms*ntypes      vel
     *
     */
    int ind=0;
    for(int i=0; i<_moltypes; ++i) {
        fread(&rec, sizeof(rec), 1, fl);  // begin record
        
        for(int a=0; a<_nmols[i]*_natoms[i]; a++) {
            float f;
            fread(&f, sizeof(float), 1, fl);
            conf.Pos(a+ind).x() = f*0.1;            
        }
        for(int a=0; a<_nmols[i]*_natoms[i]; a++) {
            float f;
            fread(&f, sizeof(float), 1, fl);
            conf.Pos(a+ind).y() = f*0.1;            
        }
        for(int a=0; a<_nmols[i]*_natoms[i]; a++) {
            float f;
            fread(&f, sizeof(float), 1, fl);
            conf.Pos(a+ind).z() = f*0.1;            
        }
        ind+=_nmols[i]*_natoms[i];
        
        I4 = rec + 1;
        while(I4!=rec)  // skip everything to end of record
            fread(&I4, sizeof(I4), 1, fl);  // end record

        if(iVel == 1) {
            fread(&rec, sizeof(rec), 1, fl);  // begin record
            
            for(int a=0; a<_nmols[i]*_natoms[i]; a++) {
                float f;
                fread(&f, sizeof(float), 1, fl);
                conf.V(a+ind).x() = f*0.1;            
            }
            for(int a=0; a<_nmols[i]*_natoms[i]; a++) {
                float f;
                fread(&f, sizeof(float), 1, fl);
                conf.V(a+ind).y() = f*0.1;            
            }
            for(int a=0; a<_nmols[i]*_natoms[i]; a++) {
                float f;
                fread(&f, sizeof(float), 1, fl);
                conf.V(a+ind).z() = f*0.1;            
            }
            
            I4 = rec + 1;
            while(I4!=rec)  // skip everything to end of record
                fread(&I4, sizeof(I4), 1, fl);  // end record
        }
    }
    if(feof(fl)) return false;
    return true;
}

