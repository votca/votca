#ifndef __VAVERAGE_H
#define	__VAVERAGE_H

#include <votca/ctp/paircalculator.h>
/**
        \brief Average charge velocity, mobility along the field, and site current.

The average velocity is calculated from occupation probabilities and charge trajectory, 
\f[ \langle \vec{v} \rangle =  \sum_{i,j}  p_j  \omega_{ji}  (\vec{r}_i - \vec{r}_j) = \hat{\mu} \vec{F}. \f]
A symmetrized local current through site \f$ i \f$ is defined as
\f[
 \vec{J}_i = \frac{1}{2} e \sum_{ j} \left( p_{j}  \omega_{ji} - p_{i} \omega_{ij} \right) \vec{r}_{ij} .
\f]
A forward local current through site \f$ i \f$ is defined as (Phys. Rev. B 79, 085203, 2009)
\f[
 \vec{J}_{i} = e \sum_{ j, r_{i, \alpha} > r_{j, \alpha} } \left( p_{j}  \omega_{ji} - p_{i} \omega_{ij} \right) \vec{r}_{ij} .
\f]
where \f$ \alpha = x,y,z \f$.

Callname: vaverage

*/

class Vaverage : public PairCalculator
{
public:
    Vaverage() {}
    ~Vaverage() {}

    const char *Description() { return "Average charge velocity, mobility along the field, and site current"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

protected:
    bool _do_mobility;
    vec _E;
    vec _v;
    int _nframes;
    map <QMCrgUnit *, vec> _symmetrized_current;
    map <QMCrgUnit *, vec> _forward_current;

};

inline void Vaverage::Initialize(QMTopology *top, Property *options) {
    _do_mobility = false;
    _nframes = 0;

    if (options->exists("options.vaverage.mobility"))
        _do_mobility = options->get("options.vaverage.mobility").as<bool > ();

    if (_do_mobility)
        _E = options->get("options.calc_rates.e_field").as<vec > ();
    _v = vec(0,0,0);
}

bool Vaverage::EvaluateFrame(QMTopology *top)
{
    PairCalculator::EvaluateFrame(top);
    _nframes++;
    return true;
}

/// Evaluate average velocity and microscopic current for each site
inline void Vaverage::EvaluatePair(QMTopology *top, QMPair *pair)
{
    QMCrgUnit *crg1 = pair->first;
    QMCrgUnit *crg2 = pair->second;

    double p1 = crg1->getOccupationProbability();
    double p2 = crg2->getOccupationProbability();

    double w12 = pair->rate12();
    double w21 = pair->rate21();
    vec r = unit<nm,m>::to(pair->r());

    _v+=(w12*p1 - w21*p2)*r;
    //cout << crg1->getOccupationProbability() << "\t" << crg2->getOccupationProbability() << "\t" << pair->rate12() <<  "\t" << pair->rate21() << "\t" <<(pair->rate12()*crg1->getOccupationProbability() - pair->rate21()*crg2->getOccupationProbability())*pair->r() << endl;

    _symmetrized_current[crg1] += 0.5 * (p2 * w21 - p1 * w12) * r;
    _symmetrized_current[crg2] += 0.5 * (p2 * w21 - p1 * w12) * r;

    if (r.getX() > 0.0) {
        _forward_current[crg1].setX(_forward_current[crg1].getX() + (p2 * w21 + p1 * w12) * r.getX());
    }
    else {
        _forward_current[crg2].setX(_forward_current[crg2].getX() - (p2 * w21 + p1 * w12) * r.getX());
    }

    if (r.getY() > 0.0) {
        _forward_current[crg1].setY(_forward_current[crg1].getY() + (p2 * w21 + p1 * w12) * r.getY());
    }
    else {
        _forward_current[crg2].setY(_forward_current[crg2].getY() - (p2 * w21 + p1 * w12) * r.getY());
    }

    if (r.getZ() > 0.0) {
        _forward_current[crg1].setZ(_forward_current[crg1].getZ() + (p2 * w21 + p1 * w12) * r.getZ());
    }
    else {
        _forward_current[crg2].setZ(_forward_current[crg2].getZ() - (p2 * w21 + p1 * w12) * r.getZ());
    }
}

/// normalize the velocity and current by the number of frames
void Vaverage::EndEvaluate(QMTopology *top)
{
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;

    for (itl = lcharges.begin(); itl != lcharges.end(); ++itl) {
      cout <<  (*itl)->getId() << "\t" 
           << (*itl)->GetCom()<< "\t" 
	   << _symmetrized_current[(*itl)] / _nframes << "\t" 
	   << _forward_current[(*itl)] / _nframes << endl;
    }

    cout << "Average velocity [m/s]: " << _v / _nframes << endl;
    if(_do_mobility)
       cout << "Average mobility [cm^2/Vs]: " << _v*_E/(_E*_E)*1E4 / _nframes << endl;
    
}

#endif	/* VAVERAGE_H */

