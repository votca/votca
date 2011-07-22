#ifndef __VAVERAGE_H
#define	__VAVERAGE_H

#include <votca/ctp/paircalculator.h>
/**
        \brief Average velocity vector and mobility (along the field) of a single charge carrier.

The average velocity is calculated from occupation probabilities and charge trajectory, 
\f[ \langle \vec{v} \rangle =  \sum_{i,j}  p_j  \omega_{ji}  (\vec{r}_i - \vec{r}_j) = \hat{\mu} \vec{F} \f].

Callname: vaverage

*/

class Vaverage : public PairCalculator
{
public:
    Vaverage() {}
    ~Vaverage() {}

    const char *Description() { return "Average velocity vector and mobility (along the field) of a single charge carrier"; }

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

protected:
    bool _do_mobility;
    vec _E;
    vec _v;
};

inline void Vaverage::Initialize(QMTopology *top, Property *options) {
    _do_mobility = false;
    if (options->exists("options.vaverage.mobility"))
        _do_mobility = options->get("options.vaverage.mobility").as<bool > ();

    if (_do_mobility)
        _E = options->get("options.calc_rates.e_field").as<vec > ();
    _v = vec(0,0,0);
}

inline void Vaverage::EvaluatePair(QMTopology *top, QMPair *pair)
{
    QMCrgUnit *crg1 = pair->first;
    QMCrgUnit *crg2 = pair->second;
    _v+=(pair->rate12()*crg1->getOccupationProbability() - pair->rate21()*crg2->getOccupationProbability())*unit<nm,m>::to(pair->r());
    //cout << crg1->getOccupationProbability() << "\t" << crg2->getOccupationProbability() << "\t" << pair->rate12() <<  "\t" << pair->rate21() << "\t" <<(pair->rate12()*crg1->getOccupationProbability() - pair->rate21()*crg2->getOccupationProbability())*pair->r() << endl;
}

void Vaverage::EndEvaluate(QMTopology *top)
{
    cout << "Average velocity [m/s]: " << _v << endl;
    if(_do_mobility)
       cout << "Average mobility [cm^2/Vs]: " << _v*_E/(_E*_E)*1E4 << endl;
}

#endif	/* VAVERAGE_H */

