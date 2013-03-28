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

#include <stdio.h>
#include "mdptopologyreader.h"
#include <vector>
#include <boost/lexical_cast.hpp>

bool MDPTopologyReader::ReadTopology(string file, Topology &top)
{
    FILE *fl;
    long I4;
    float R4;
    double R8;
        
    fl = fopen(file.substr(0, file.length()-4).c_str(), "rb");
    if(!fl)
        return false;
    
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
    cout << "I have " << (I4 ? "" : "no ") << "velocities\n";
    
    double dt, unit;
    vec box;
    long moltypes;
    
    fread(&dt, sizeof(dt), 1, fl);
    fread(&box, sizeof(box), 1, fl);
    fread(&unit, sizeof(unit), 1, fl);
    fread(&moltypes, sizeof(moltypes), 1, fl);
    
    cout << dt << " " << box << " " << unit << " " << moltypes << endl;
    
    vector<long> nmols;
    vector<long> natoms;
    
    for(int i=0; i<moltypes; ++i) {
        fread(&I4, sizeof(I4), 1, fl);
        printf("%X\n", I4);
        nmols.push_back(I4);
        cout << "nmols " << I4 << endl;
        fread(&I4, 1, sizeof(I4), fl);
        natoms.push_back(I4);
        cout << "natoms " << I4 << endl;
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
      
     
    for(int type=0; type<moltypes; ++type) {
        vector<double> masses;             
        for(int atom=0; atom<natoms[type]; ++atom) {
            fread(&R8, sizeof(R8), 1, fl);
            masses.push_back(R8);
        }
        // is molecule read in trajectory?
        fread(&I4, sizeof(I4), 1, fl);
        // we create it anyway...
        for(int mol=0; mol<nmols[type]; ++mol) {
            MoleculeInfo *mi = top.CreateMolecule(boost::lexical_cast<string>(type+1));
            for(int atom=0; atom<natoms[type]; ++atom) {
                mi->AddBead(top.CreateBead(1, "", top.GetOrCreateBeadType("no"), 0, R8, 0)->getId(),
                  boost::lexical_cast<string>(atom+1));                                
            }
        }
    }
    
    // read record begin
    fread(&I4, sizeof(I4), 1, fl);
    
    fclose(fl);
    
    return true;
}
