#include <votca/ctp/polarsite.h>

namespace votca { namespace ctp {


void PolarSite::ImportFrom(PolarSite *templ) {

    _pos = templ->getPos();
    


}


}}

