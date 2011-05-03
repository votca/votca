#ifndef __VOTCA_MD2QM_AVGVELOCITY_H
#define	__VOTCA_MD2QM_AVGVELOCITY_H

#include <paircalculator.h>
class AvgVelocity : public PairCalculator
{
public:
    AvgVelocity() {}
    ~AvgVelocity() {}

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
    _v+=(pair->rate12()*crg1->getOccupationProbability() - pair->rate21()*crg1->getOccupationProbability())*pair->r();
}

void AvgVelocity::EndEvaluate(QMTopology *top)
{
    cout << "Average velocity: " << _v << endl;
    if(_do_mobility)
       cout << "Average mobility: " << _v*_E/(_E*_E) << endl;
}

#endif	/* AVGVELOCITY_H */

