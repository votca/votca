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

#ifndef _cgengine_H
#define	_cgengine_H

#include <list>
#include <map>
#include <boost/program_options.hpp>
#include "topology.h"
#include "cgmoleculedef.h"
#include <votca/tools/datacollection.h>
#include "topologymap.h"
#include "cgobserver.h"

#include <votca/tools/vec.h>
#include "cgmoleculedef.h"
#include "cgengine.h"
#include "molecule.h"
#include "topologyreader.h"
#include "trajectorywriter.h"
#include "trajectoryreader.h"
#include <votca/tools/tokenizer.h>
#include <votca/tools/matrix.h>
#include "nematicorder.h"

using namespace std;

/**
    \brief coarse graining engine

    This class manages the coarse graining, at the moment it does the measurement stuff

    TODO: split this into an additional VotcaApplication object

*/
class CGEngine
{
public:
    CGEngine();
    ~CGEngine();    

    void Initialize();
    void ParseCommandLine(int argc, char **argv);

    /**
        create a coarse grained topolgy based on a given topology
    */
    TopologyMap *CreateCGTopology(Topology &in, Topology &out);
    
    /**
        load molecule type from file
    */
    void LoadMoleculeType(string filename);
    
    
    /**
        begin the coarse graining process
     */
    void BeginCG(Topology &top);
    /**
        ends the coarse graining process
     */
    void EndCG();
          
    /**
        evaluate configuration while in coarse graining process
     */
    void EvalConfiguration(Topology &top_cg);
        
    CGMoleculeDef *getMoleculeDef(string name);
    
    void AddObserver(CGObserver *observer);

    boost::program_options::options_description_easy_init
        AddProgramOptions() { return _op_desc_specific.add_options(); }

    boost::program_options::variables_map &OptionsMap() { return _op_vm; }
    boost::program_options::options_description &OptionsDesc() { return _op_desc; }


    void Run();
    
private:
    void InitializeStandardOptions();

    list<CGObserver *> _observers;
    map<string, CGMoleculeDef *> _molecule_defs;

    boost::program_options::options_description _op_desc;
    boost::program_options::options_description _op_desc_specific;
    boost::program_options::variables_map _op_vm;
};

inline CGMoleculeDef *CGEngine::getMoleculeDef(string name)
{
    map<string, CGMoleculeDef*>::iterator iter;
    
    // if there is only 1 molecule definition, don't care about the name
    if(_molecule_defs.size()==1 && name == "unnamed") {
            return (*(_molecule_defs.begin())).second;
    }
    
    iter = _molecule_defs.find(name);        
    if(iter == _molecule_defs.end()) return NULL;
    return (*iter).second;
}

inline void CGEngine::AddObserver(CGObserver *observer)
{
    _observers.push_back(observer);
}

#endif	/* _cgengine_H */

