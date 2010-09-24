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

#include "lammpsreader.h"

namespace votca { namespace csg {

using namespace std;

bool LAMMPSReader::ReadTopology(string file,  Topology &top)
{
    ifstream fl;
    fl.open(file.c_str());
    if(!fl.is_open())
        throw std::ios_base::failure("Error on open topologyl file: " + file);

    fl.close();

}

bool LAMMPSReader::Open(const string &file)
{
    _fl.open(file.c_str());
    if(!_fl.is_open())
        throw std::ios_base::failure("Error on open topologyl file: " + file);
}

void LAMMPSReader::Close()
{
    _fl.close();
}

bool LAMMPSReader::FirstFrame(Topology &conf)
{
    return true;
}

bool LAMMPSReader::NextFrame(Topology &conf)
{
    return true;
}

}}
