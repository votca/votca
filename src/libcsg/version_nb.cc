#include <votca/tools/version.h>
#include <iostream>
#include "version.h"

namespace votca { namespace csg {

static const std::string version_str = "VERSION NOT SET (compiled " __DATE__ ", " __TIME__ ")";


const std::string &CsgVersionStr()
{
    return version_str;
}

void HelpTextHeader(const std::string &tool_name)
{
    std::cout <<
         << "\tVOTCA ( http://www.votca.org )\n"
         << tool_name << ", version " << votca::csg::CsgVersionStr()
         << "\nvotca_tools, version " << votca::tools::ToolsVersionStr() 
         << "\n\n";
}

}}

