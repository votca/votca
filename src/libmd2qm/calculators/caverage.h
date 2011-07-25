#ifndef CAVERAGE_H
#define CAVERAGE_H

#include <votca/ctp/qmpair.h>
#include <votca/ctp/paircalculator.h>

/**                                  
    \brief Average microscopic current.         

Callname: caverage 

Average value of the microscopic current. A local current through site $i$ is defined as
\f[
 \vec{J}_i = \frac{1}{2} e \sum_{ j} \left( p_{j}  \omega_{ji} - p_{i} \omega_{ij} \right) \vec{r}_{ij} .
\f]

*/                                         

class Caverage : public PairCalculator
{
public:
    Caverage() {};
    ~Caverage() {};

    const char *Description() { return "Average microscopic current."; }

    void EvaluatePair(QMTopology *top, QMPair *pair);
    void Initialize(QMTopology *top, Property *options);
    void EndEvaluate(QMTopology *top);
    bool EvaluateFrame(QMTopology *top);

private:
    //frame counter
    int _nframes;
    map <QMCrgUnit *, vec> _map_crgunt_current;
    map <QMCrgUnit *, vec> _map_crgunt_current2;
};


void Caverage::Initialize(QMTopology *top, Property *options){
    _nframes = 0;
    
}


bool Caverage::EvaluateFrame(QMTopology *top)
{
    // evaluate frame
    PairCalculator::EvaluateFrame(top);
    _nframes++;
    return true;
}

inline void Caverage::EvaluatePair(QMTopology *top, QMPair *pair){
    QMCrgUnit *crg1 = pair->Crg1();
    QMCrgUnit *crg2 = pair->Crg2();

    double p1 = crg1->getOccupationProbability();
    double p2 = crg2->getOccupationProbability();

    double w12 = pair->rate12();
    double w21 = pair->rate21();
    vec r = pair->r();

    _map_crgunt_current[crg1] += 0.5 * (p2 * w21 - p1 * w12) * r;
    _map_crgunt_current[crg2] += 0.5 * (p2 * w21 - p1 * w12) * r;

    if (r.getX() > 0.0) {
        _map_crgunt_current2[crg1].setX(_map_crgunt_current2[crg1].getX() + (p2 * w21 + p1 * w12) * r.getX());
    }
    else {
        _map_crgunt_current2[crg2].setX(_map_crgunt_current2[crg2].getX() - (p2 * w21 + p1 * w12) * r.getX());
    }

    if (r.getY() > 0.0) {
        _map_crgunt_current2[crg1].setY(_map_crgunt_current2[crg1].getY() + (p2 * w21 + p1 * w12) * r.getY());
    }
    else {
        _map_crgunt_current2[crg2].setY(_map_crgunt_current2[crg2].getY() - (p2 * w21 + p1 * w12) * r.getY());
    }

    if (r.getZ() > 0.0) {
        _map_crgunt_current2[crg1].setZ(_map_crgunt_current2[crg1].getZ() + (p2 * w21 + p1 * w12) * r.getZ());
    }
    else {
        _map_crgunt_current2[crg2].setZ(_map_crgunt_current2[crg2].getZ() - (p2 * w21 + p1 * w12) * r.getZ());
    }


}


void Caverage::EndEvaluate(QMTopology* top) {
    // normalize
    vector<QMCrgUnit *> lcharges = top->CrgUnits();
    vector<QMCrgUnit *>::iterator itl;
    for (itl = lcharges.begin(); itl!=lcharges.end(); ++itl) {
      cout <<  (*itl)->getId() << "\t" << (*itl)->GetCom()<< "\t" << _map_crgunt_current[(*itl)] / _nframes <<
              "\t" << _map_crgunt_current2[(*itl)] / _nframes << endl;
    }
}

#endif	/* CAVERAGE_H */

