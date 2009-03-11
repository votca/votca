/* 
 * File:   pair.h
 * Author: ruehle
 *
 * Created on March 5, 2009, 10:30 AM
 */

#ifndef _BEADPAIR_H
#define	_BEADPAIR_H

/**
   \brief A particle pair
 
   This class defines a particle pair. The future plan is, that the Pair class
   can be overloaded and Particle list creates these inherited pairs.
 
 */

class BeadPair
    : public std::pair<Bead *, Bead *>
{
public:
    BeadPair() {}
    BeadPair(Bead *bead1, Bead *bead2, vec r)
            : std::pair<Bead *, Bead *>(bead1, bead2), _r(r), _dist(abs(r)) {}
        
    vec &r() { return _r; }     
    double &dist() { return _dist; }
    
private:
        vec _r;
        double _dist;
};

#endif	/* _PAIR_H */

