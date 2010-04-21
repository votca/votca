/* 
 * File:   ctp_observer.h
 * Author: vehoff
 *
 * Created on April 8, 2010, 2:09 PM
 */

#include "qmtopology.h"

#ifndef _CTP_OBSERVER_H
#define	_CTP_OBSERVER_H

class CTPObserver
{
public:
    /// \brief called before the first frame
    virtual void BeginCTP(QMTopology *top) = 0;
    /// \brief called after the last frame
    virtual void EndCTP(QMTopoplogy *top) = 0;
    // \brief called for each frame
    virtual void EvalCTP(QMTopology *top) = 0;
};

#endif	/* _CTP_OBSERVER_H */

