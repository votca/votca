/// \addtogroup csg
///@{
// 
// File:   gmxtrajectoryreader.h
// Author: ruehle
//
// Created on April 5, 2007, 2:41 PM
//

#ifndef _gmxtrajectoryreader_H
#define	_gmxtrajectoryreader_H

#include <string>
#include "trajectoryreader.h"

using namespace std;

namespace gmx {
   extern "C" {
        #include <statutil.h>
        #include <typedefs.h>
        #include <smalloc.h>
        #include <vec.h>
        #include <copyrite.h>
        #include <statutil.h>
        #include <tpxio.h>
    }
    // this one is needed because of bool is defined in one of the headers included by gmx
    #undef bool
}

/**
    \brief class for reading gromacs trajectory files

    This class provides the TrajectoryReader interface and encapsulates the trajectory reading function of gromacs

*/
class GMXTrajectoryReader : public TrajectoryReader
{
    public:        
        /// open a trejectory file
        bool Open(const string &file);
        /// read in the first frame
        bool FirstFrame(Configuration &conf);
        /// read in the next frame
        bool NextFrame(Configuration &conf);

        void Close();
        
        TrajectoryReader *Clone() { return dynamic_cast<TrajectoryReader*>(new GMXTrajectoryReader()); }

    private:
        string _filename;
        
        // gmx status used in read_first_frame and _read_next_frame;
        int _gmx_status;
        /// gmx frame
        gmx::t_trxframe _gmx_frame;
        
};

#endif	/* _gmxtrajectoryreader_H */

/// @}
