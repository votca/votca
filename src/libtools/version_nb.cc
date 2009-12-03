// This file is for netbeans only
#include "version.h"

namespace votca { namespace tools {

static const std::string version_str = "VERSION NOT SET (compiled " __DATE__ ", " __TIME__ ")";

const std::string &ToolsVersionStr()
{
    return version_str;
}

}}

