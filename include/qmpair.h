/* 
 * File:   qmbeadpair.h
 * Author: james
 *
 * Created on 14 December 2009, 08:47
 */

#ifndef _QMPAIR_H
#define	_QMPAIR_H

#include <moo/crgunit.h>
#include <utility>

class QMTopology;

class QMPair :
    public std::pair<CrgUnit *, CrgUnit *>
{
public:
    QMPair() {}
    QMPair(CrgUnit *crg1, CrgUnit *crg2, QMTopology * top);

    ~QMPair(){
        if(_ghost != NULL)
            delete _ghost;
    }
    /// \brief the vector connecting two beads
    vec &r() { return _r; }
    /// \brief the distance of the beads
    double &dist() { return _dist; }
    /// \brief returns the transfer integral
    double &j(){ return _j;}
    /// \brief returns the transfer rate
    double &rate(){return _rate;}
    /// \brief set the transfer integral
    void setJ(const double & j){_j=j;}
    /// \brief set the transfer integral as the rms of the transfer integrals
    void setJ(vector <double> js);
    /// \brief set the transfer rate
    void setRate(double rate) {_rate=rate;}

protected:
    /// vector connecting the two beads
    vec _r;
    /// distance between the two beads
    double _dist;
    /// transfer integral
    double _j;
    /// transfer rate
    double _rate;
    /// ghost atom in case the molecules are neighbors across a boundary
    CrgUnit * _ghost;
};

#endif	/* _QMBEADPAIR_H */

