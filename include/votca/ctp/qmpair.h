/* 
 * File:   qmbeadpair.h
 * Author: james
 *
 * Created on 14 December 2009, 08:47
 */

#ifndef _QMPAIR_H
#define	_QMPAIR_H

#include "qmcrgunit.h"
#include "customfields.h"
#include <utility>

namespace votca { namespace ctp {

class QMTopology;

class QMPair :
    public std::pair<QMCrgUnit *, QMCrgUnit *>, public CustomFields
{
public:
    QMPair(): _r(0.,0.,0.), _rate_12(0.),_rate_21(0.),_ghost(NULL),_crg2(NULL), _in_database(false) {}

    QMPair(QMCrgUnit *crg1, QMCrgUnit *crg2, QMTopology * top);

    ~QMPair(){
        if(_ghost != NULL)
            delete _ghost;
    }
    /**
     * \brief the vector connecting two beads
     * @return pbc correct distance vector
     */
    vec &r() { return _r; }

    /**
     * \brief the distance of the beads
     * @return pbc correct distance
     */
    double dist() { return abs(_r); }

    /**
     * \brief vector of transfer integrals (vector because of degenerate orbitals)
     * @return vector of transfer integrals
     */
    vector <double> &Js() { return _Js; }

    /**
     * \brief transfer rate from first to second
     * @return transfer rate from first to second
     */

    double &rate12(){return _rate_12;}

    /**
     * \brief transfer rate from second to first
     * @return transfer rate from second to first
     */     
    double &rate21(){return _rate_21;}

    /**
     * \brief set the transfer integral
     * @param js vector with transfer integrals
     */
    void setJs(const vector <double> &js){_Js=js;}

    /** \brief calculate the effective transfer integral as the rms of the transfer integrals
     * @return effective tranfer integral squared
     */
    double calcJeff2();

    /**
     *  \brief set transfer rate from first to second
     * @param rate  forward rate
     */
    void setRate12(double rate) {_rate_12=rate;}
    
    /**
     * \brief set transfer rate from second to first
     * @param r backward rate
     */
    void setRate21(double rate) {_rate_21=rate;}

    /**
     * \brief first crg unit (might be ghost copy for pbc image)
     *
     * first and second are original crg units, crg1 and crg2 take into account
     * pbc and might be ghost copies. If you don't need the PBC corrected
     * positions, never (!!) use this function, especially if properties of
     * a charge unit are changed. Use pair.first instead!

     * @return crg unit 1
     */
    QMCrgUnit *Crg1PBCCopy() {return first;}

    /**
     * \brief second crg unit (might be ghost copy for pbc image)
     *
     * first and second are original crg units, crg1 and crg2 take into account
     * pbc and might be ghost copies. If you don't need the PBC corrected
     * positions, never (!!) use this function, especially if properties of
     * a charge unit are changed. Use pair.second instead!
     *
     * @return crg unit 2
     */
    QMCrgUnit *Crg2PBCCopy() {return _crg2;}

    void setInDatabase(bool indb) { _in_database = indb; }
    bool getInDatabase() { return _in_database; }

    int getId() { return _id; }
    void setId(int id) { _id = id; }

protected:
    /// vector connecting the two beads
    vec _r;
    /// transfer integrals, multiple entries in case of degeneracy
    vector <double> _Js;
    /// transfer rate from first to second
    double _rate_12;
    /// transfer rate from second to first
    double _rate_21;
    
    /// ghost atom in case the molecules are neighbors across a boundary
    QMCrgUnit * _ghost;

    QMCrgUnit *_crg2;

    bool _in_database;
    int _id;
};

}}

#endif	/* _QMBEADPAIR_H */

