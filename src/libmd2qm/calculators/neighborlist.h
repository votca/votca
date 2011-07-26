#ifndef __NEIGHBORLIST_H
#define	__NEIGHBORLIST_H

#include <votca/ctp/qmcalculator.h>

/** \brief Constructs a list of neighboring conjugated segments.

Callname: neighborlist

Two segments are added to this list if the distance between centers of mass of any their rigid fragments is below a certain cutoﬀ. This allows neighbors to be selected on a criterion of minimum distance of approach rather than center of mass distance, which is useful for molecules with anisotropic shapes.

*/

class Neighborlist : public QMCalculator
{
public:
    Neighborlist() {}
    ~Neighborlist() {}

    const char *Description() { return "Constructs a list of neighboring conjugated segments"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    
protected:
    double _cutoff;
};

inline void Neighborlist::Initialize(QMTopology *top, Property *options){
    _cutoff = options->get("options.neighborlist.cutoff").as<double>();
}

inline bool Neighborlist::EvaluateFrame(QMTopology *top)
{
    top->nblist().Cleanup();
    top->nblist().setCutoff(_cutoff);
    BeadList list;
    list.Generate(*top, "*");
    top->nblist().Generate(list);
}


#endif	/* __NEIGHBORLIST_H */

