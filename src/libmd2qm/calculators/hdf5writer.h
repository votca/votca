#ifndef __VOTCA_MD2QM_HDF5WRITE_H
#define __VOTCA_MD2QM_HDF5WRITE_H

#include "qmcalculator.h"
#include "statesaverh5.h"

class HDF5Writer : public QMCalculator {
public:
    PairCalculator() {};
    virtual ~PairCalculator() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);

private:
    StateSaverH5 _saver;
};

inline void HDF5Writer::Initialize(QMTopology *top, Property *options)
{
    _saver.Open(*top, "state.h5");
}

inline bool HDF5Writer::EvaluateFrame(QMTopology *top)
{
    _saver.WriteFrame();
    return true;
}

inline void HDF5Writer::EndEvaluate(QMTopology *top)
{
    _saver.Close();
}

#endif

