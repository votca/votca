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

#ifndef __VOTCA_CTP_TOPOLOGY_H
#define	__VOTCA_CTP_TOPOLOGY_H

#include <votca/tools/property.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/fragment.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/molecule.h>

namespace votca { namespace ctp {

/**
    \brief Container for all types of conjugated segments and rigid fragments

    Given a coarse-grained topology, it re-orients the and shifts rigid fragments. 
*/
class Topology 
{
public:
    Topology();
   ~Topology();

    /// Load the topology based on definitions of conjugated segments
    void Initialize( Property &topology );

    int getDatabaseId() { return _db_id; };
    void setDatabaseId(int id) { _db_id = id; }

protected:

    vector < Atom* > _atoms;
    vector < Fragment* > _fragments;
    vector < Segment* > _segments;
    vector < Molecule* > _molecules;
    
    map < int, Segment* > _map_id_segment;
    map < string, Segment* > _map_name_segment;

    int _db_id;
};

}}

#endif	/* __VOTCA_CTP_TOPOLOGY_H */

