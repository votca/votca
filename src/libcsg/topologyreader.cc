// 
// File:   topologyreader.cc
// Author: ruehle
//
// Created on January 18, 2008, 6:12 PM
//

#include "topologyreader.h"
#include "modules/io/gmxtopologyreader.h"
#include "modules/io/grotopologyreader.h"
#include "modules/io/xmltopologyreader.h"
#include "modules/io/pdbtopologyreader.h"

void TopologyReader::RegisterPlugins(void)
{
    TopReaderFactory().Register("tpr", new GMXTopologyReader(), false);
    TopReaderFactory().Register("gro", new GROTopologyReader(), false); 
    TopReaderFactory().Register("xml", new XMLTopologyReader(), false); 
    TopReaderFactory().Register("pdb", new PDBTopologyReader(), false); 
//    TopReaderFactory().Register("mdp", new MDPTopologyReader(), false); 
}
