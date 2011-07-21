#ifndef __VOTCA_MD2QM_AVGVELOCITY_H
#define	__VOTCA_MD2QM_AVGVELOCITY_H

#include <votca/ctp/paircalculator.h>

class AvgVelocity : public PairCalculator
{
public:
    AvgVelocity() {}
    ~AvgVelocity() {}

    const char *Description() { return "TODO"; }

    void Initialize(QMTopology *top, Property *options);
    void EvaluatePair(QMTopology *top, QMPair *pair);
    void EndEvaluate(QMTopology *top);

protected:
    bool _do_mobility;
    vec _E;
    vec _v;
};

inline void AvgVelocity::Initialize(QMTopology *top, Property *options) {
    _do_mobility = false;
    if (options->exists("options.avgvelocity.mobility"))
        _do_mobility = options->get("options.avgvelocity.mobility").as<bool > ();

    if (_do_mobility)
        _E = options->get("options.calc_rates.e_field").as<vec > ();
    _v = vec(0,0,0);
}

inline void AvgVelocity::EvaluatePair(QMTopology *top, QMPair *pair)
{
    QMCrgUnit *crg1 = pair->first;
    QMCrgUnit *crg2 = pair->second;
    _v+=(pair->rate12()*crg1->getOccupationProbability() - pair->rate21()*crg2->getOccupationProbability())*unit<nm,m>::to(pair->r());
    //cout << crg1->getOccupationProbability() << "\t" << crg2->getOccupationProbability() << "\t" << pair->rate12() <<  "\t" << pair->rate21() << "\t" <<(pair->rate12()*crg1->getOccupationProbability() - pair->rate21()*crg2->getOccupationProbability())*pair->r() << endl;
}

void AvgVelocity::EndEvaluate(QMTopology *top)
{
    cout << "Average velocity [m/s]: " << _v << endl;
    if(_do_mobility)
       cout << "Average mobility [cm^2/Vs]: " << _v*_E/(_E*_E)*1E4 << endl;
}

#endif	/* AVGVELOCITY_H */

