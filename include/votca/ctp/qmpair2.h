#ifndef _QMPAIR2_H
#define _QMPAIR2_H

#include "segment.h"
#include <utility>


namespace votca { namespace ctp {

class Topology;



class QMPair2 : public std::pair< Segment*, Segment* >
{
public:
    QMPair2() { };
   ~QMPair2() { };

private:

    int _id;
};





}}




#endif
