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
    TopReaderFactory().Register<GMXTopologyReader>("tpr");
    TopReaderFactory().Register<GROTopologyReader>("gro");
    TopReaderFactory().Register<XMLTopologyReader>("xml");
    TopReaderFactory().Register<PDBTopologyReader>("pdb");
}
