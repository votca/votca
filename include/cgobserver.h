// 
// File:   cgobserver.h
// Author: victor
//
// Created on 4. Juni 2008, 17:08
//

#ifndef _CGOBSERVER_H
#define	_CGOBSERVER_H

#include "topology.h"

/**
   \brief Observer class for analysis hook

   Each application which performs analysis operations should use CGEngine. It offers
   a hook (callback class) during the coarse-graining process to evaluate each frame.
   The user does not have to take care about mapping and other stoff. Just oberload
   this class and analyze properties of interest.

 */

class CGObserver
{
public:
    /// \brief called before the first frame
    virtual void BeginCG(Topology *top, Topology *top_atom = 0) = 0;
    /// \brief called after the last frame
    virtual void EndCG() = 0;
    // \brief called for each frame which is mapped
    virtual void EvalConfiguration(Topology *top, Topology *top_atom = 0) = 0;
};


#endif	/* _CGOBSERVER_H */

