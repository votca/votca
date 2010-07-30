/* 
 * File:   energycorr.h
 * Author: lukyanov
 *
 * Created on July 28, 2010, 6:04 PM
 */

#ifndef ENERGYCORR_H
#define	ENERGYCORR_H

#include "qmpair.h"
#include "qmcalculator.h"
#include <votca/tools/histogramnew.h>

/// Calculate reduced site energy correlation function, see J.Chem.Phys. 129, 034709 (2008)

class EnergyCorr : public QMCalculator
{
public:
    EnergyCorr() {};
    ~EnergyCorr() {};

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void EndEvaluate(QMTopology *top);
    
private:
    // histogram for recording energy correlation function
    HistogramNew _hist_corr;
    // histogram for recording counts (needed to normalize the first one)
    HistogramNew _hist_count;
    // histogram for recording correlation function squared (needed to calculate error bars)
    HistogramNew _hist_corr2;
    double _min, _max;
    int _nbins;
    /// name of the output file
    string _outfile;


};

#endif	/* ENERGYCORR_H */

