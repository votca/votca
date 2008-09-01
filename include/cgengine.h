/// \addtogroup csg
///@{
// 
// File:   cgengine.h
// Author: ruehle
//
// Created on April 23, 2007, 6:01 PM
//

#ifndef _cgengine_H
#define	_cgengine_H

#include <list>
#include <map>
#include <boost/program_options.hpp>
#include "topology.h"
#include "configuration.h"
#include "cgmoleculedef.h"
#include <tools/datacollection.h>
#include "topologymap.h"
#include "cgobserver.h"

#include <tools/vec.h>
#include "cgmoleculedef.h"
#include "cgengine.h"
#include "molecule.h"
#include "topologyreader.h"
#include "trajectorywriter.h"
#include "trajectoryreader.h"
#include <tools/tokenizer.h>
#include <tools/matrix.h>
#include "nematicorder.h"

using namespace std;

/**
    \brief coarse graining engine

    This class manages the coarse graining, at the moment it does the measurement stuff

*/
class CGEngine
{
public:
    CGEngine() {}
    ~CGEngine();    
    
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
    void EvalConfiguration(Configuration &conf_cg);
        
    CGMoleculeDef *getMoleculeDef(string name);
    
    void AddObserver(CGObserver *observer);
    
    
    void AddProgramOptions(boost::program_options::options_description &desc);
    void Run(boost::program_options::options_description &desc, boost::program_options::variables_map &vm);
    
private:
    list<CGObserver *> _observers;
    map<string, CGMoleculeDef *> _molecule_defs;
};

inline CGMoleculeDef *CGEngine::getMoleculeDef(string name)
{
    map<string, CGMoleculeDef*>::iterator iter;
    
    // if there is only 1 molecule definition, don't care about the name
    if(_molecule_defs.size()==1) {
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

/// @}
