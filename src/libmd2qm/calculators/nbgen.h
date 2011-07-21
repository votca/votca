#ifndef __VOTCA_MD2QM_NBGEN_H
#define	__VOTCA_MD2QM_NBGEN_H

#include <votca/ctp/qmcalculator.h>

class NBGen : public QMCalculator
{
public:
    NBGen() {}
    ~NBGen() {}

    const char *Description() { return "Regenerate the neigbour list with a new cutoff"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    
protected:
    double _cutoff;
};

inline void NBGen::Initialize(QMTopology *top, Property *options){
    _cutoff = options->get("options.nbgen.cutoff").as<double>();
}

inline bool NBGen::EvaluateFrame(QMTopology *top)
{
    top->nblist().Cleanup();
    top->nblist().setCutoff(_cutoff);
    BeadList list;
    list.Generate(*top, "*");
    top->nblist().Generate(list);


}


#endif	/* NBGEN_H */

