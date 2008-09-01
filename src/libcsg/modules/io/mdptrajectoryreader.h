/// \addtogroup csg
///@{
/* 
 * File:   mdptrajectoryreader.h
 * Author: victorr
 *
 * Created on June 12, 2008, 1:51 PM
 */

#ifndef _MDPTRAJECTORYREADER_H
#define	_MDPTRAJECTORYREADER_H

#include <stdio.h>
#include "trajectoryreader.h"

class MDPTrajectoryReader : public TrajectoryReader
{
    public:        
        /// open a trejectory file
        bool Open(const string &file);
        /// read in the first frame
        bool FirstFrame(Configuration &conf);
        /// read in the next frame
        bool NextFrame(Configuration &conf);

        void Close();
        
        TrajectoryReader *Clone() { return dynamic_cast<TrajectoryReader*>(new MDPTrajectoryReader()); }

    private:
        FILE *_fl;
        int _moltypes;
        vector<int> _nmols;
        vector<int> _natoms;
};



#endif	/* _MDPTRAJECTORYREADER_H */

/// @}
